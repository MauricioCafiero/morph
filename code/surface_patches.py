#!/usr/bin/env python3
"""
surface_patches.py — compare charged / hydrophobic surface patches between two
conformations of the same peptide (e.g. apo fold vs Boltz co-fold).

For each state it computes, per residue: formal charge, Kyte-Doolittle
hydrophobicity, and relative solvent-accessible surface area (rSASA, mdtraj).
An exposed residue (rSASA > --expose) that is charged or strongly hydrophobic is
a "spot"; spatially adjacent spots of the same type are grouped into patches.
Then it diffs the two states: which spots persist, bury, newly expose, and which
patches are gained or lost when the peptide folds into the bound form.

Outputs (into --outdir):
  <tag>_<state>.pdb   peptide, B-factor = residue charge, occupancy = hydrophobicity
                      (PyMOL: `spectrum b, blue_white_red` for charge;
                              `spectrum q, white_orange`  for hydrophobicity)
  <tag>_residues.csv  per-residue charge / KD / rSASA(apo,cof) / exposure change
  <tag>_summary.txt   patches per state and the apo->cofold diff

Reads .pdb and Boltz/OpenFold .cif (peptide ATOM records only) via morph.py.

Usage:
  surface_patches.py --tag s3_esm2_f4 --outdir patches apo.pdb cofold.cif
  surface_patches.py --labels apo,boltz --expose 0.25 ... apo.pdb cofold.cif
"""
import argparse
import os
import sys
import numpy as np
import mdtraj as md

import morph  # same folder: parse_pdb / parse_cif / kabsch / pdb_line

# Kyte-Doolittle hydrophobicity (higher = more hydrophobic)
KD = {"ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
      "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3,
      "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5,
      "LYS": -3.9, "ARG": -4.5}
CHARGE = {"ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1}        # HIS treated neutral at pH 7
# Tien et al. 2013 theoretical max SASA per residue (A^2), for relative SASA
MAXSASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLU": 223,
           "GLN": 225, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
           "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
           "TYR": 263, "VAL": 174}
HYDRO_CUT = 1.5     # KD above this = "strongly hydrophobic"
PATCH_CUT = 8.0     # CB-CB distance (A) linking two spots into one patch


def write_peptide_pdb(atoms_by_key, keys, path):
    """Emit a peptide-only heavy-atom PDB from a morph atom dict."""
    with open(path, "w") as fh:
        for serial, k in enumerate(keys, 1):
            name, resname, resseq, element, xyz = atoms_by_key[k]
            fh.write(morph.pdb_line(serial, name, resname, resseq, element, xyz, "A"))
        fh.write("END\n")


def per_residue(pdb_path):
    """Return (resnums, resnames, rSASA[%0-1], CB_coords[A]) for a peptide PDB."""
    t = md.load(pdb_path)
    sasa = md.shrake_rupley(t, mode="residue")[0] * 100.0        # nm^2 -> A^2
    resnums, resnames, rsasa, cb = [], [], [], []
    for i, res in enumerate(t.topology.residues):
        rn = res.name
        resnums.append(res.resSeq)
        resnames.append(rn)
        rsasa.append(sasa[i] / MAXSASA.get(rn, 200.0))
        atom = next((a for a in res.atoms if a.name == "CB"), None) \
            or next((a for a in res.atoms if a.name == "CA"), None)
        cb.append(t.xyz[0, atom.index] * 10.0)                   # nm -> A
    return np.array(resnums), resnames, np.array(rsasa), np.array(cb)


def find_patches(idxs, coords, cutoff=PATCH_CUT):
    """Connected components of spot residues whose CB are within cutoff."""
    idxs = list(idxs)
    seen, patches = set(), []
    for s in idxs:
        if s in seen:
            continue
        stack, comp = [s], []
        while stack:
            a = stack.pop()
            if a in seen:
                continue
            seen.add(a)
            comp.append(a)
            for b in idxs:
                if b not in seen and np.linalg.norm(coords[a] - coords[b]) <= cutoff:
                    stack.append(b)
        patches.append(sorted(comp))
    return sorted(patches, key=len, reverse=True)


def analyze(structures, tag="peptide", labels=("apo", "cofold"), outdir="patches",
            chain="A", expose=0.25, verbose=True):
    """Compare charged/hydrophobic surface patches between two states. Returns a summary dict."""
    labels = list(labels)
    os.makedirs(outdir, exist_ok=True)

    # load both, keep atoms common to both in the first structure's order
    structs = [morph.load(p, chain) for p in structures]
    common = set(structs[0]) & set(structs[1])
    keys = [k for k in structs[0] if k in common]
    if not keys:
        sys.exit("no shared atoms between the two structures")

    # superpose cofold onto apo (CA) so the two colored structures overlay in PyMOL
    ca = np.array([k[1] == "CA" for k in keys])
    P = np.array([structs[1][k][4] for k in keys])
    Q = np.array([structs[0][k][4] for k in keys])
    R, t = morph.kabsch(P[ca], Q[ca])
    for n, k in enumerate(keys):
        name, resname, resseq, element, _ = structs[1][k]
        structs[1][k] = (name, resname, resseq, element, tuple(P[n] @ R.T + t))

    # write peptide PDBs and run per-residue analysis on each
    data = []
    for st, lab in zip(structs, labels):
        p = os.path.join(outdir, f"{tag}_{lab}.pdb")
        write_peptide_pdb(st, keys, p)
        data.append(per_residue(p))
    (rn_a, names_a, sasa_a, cb_a) = data[0]
    (rn_b, names_b, sasa_b, cb_b) = data[1]
    assert names_a == names_b, "residue sequences differ between the two states"
    n = len(names_a)
    charge = np.array([CHARGE.get(r, 0) for r in names_a])
    kd = np.array([KD.get(r, 0.0) for r in names_a])

    # bake charge -> B-factor, hydrophobicity -> occupancy, for PyMOL surface coloring
    for st, lab in zip(structs, labels):
        p = os.path.join(outdir, f"{tag}_{lab}.pdb")
        lines = []
        for serial, k in enumerate(keys, 1):
            name, resname, resseq, element, xyz = st[k]
            i = int(np.where(rn_a == resseq)[0][0]) if resseq in rn_a else 0
            nm = name if len(name) >= 4 else f" {name:<3}"
            q = max(min(kd[i] / 4.5, 1.0), 0.0)                 # 0..1 hydrophobicity -> occupancy
            lines.append(f"ATOM  {serial:5d} {nm} {resname:>3} A{resseq:4d}   "
                         f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{q:6.2f}{charge[i]:6.2f}"
                         f"          {element:>2}\n")
        open(p, "w").writelines(lines + ["END\n"])

    # classify spots per state
    def spots(sasa):
        exposed = sasa > expose
        pos = [i for i in range(n) if exposed[i] and charge[i] > 0]
        neg = [i for i in range(n) if exposed[i] and charge[i] < 0]
        pho = [i for i in range(n) if exposed[i] and kd[i] >= HYDRO_CUT]
        return exposed, pos, neg, pho
    exp_a, pos_a, neg_a, pho_a = spots(sasa_a)
    exp_b, pos_b, neg_b, pho_b = spots(sasa_b)

    def rtag(i):
        return f"{names_a[i]}{rn_a[i]}"

    # per-residue CSV
    csv = os.path.join(outdir, f"{tag}_residues.csv")
    with open(csv, "w") as fh:
        fh.write(f"resnum,resname,charge,kd,rSASA_{labels[0]},rSASA_{labels[1]},d_rSASA,class,change\n")
        for i in range(n):
            cls = "pos" if charge[i] > 0 else "neg" if charge[i] < 0 else ("hydrophobic" if kd[i] >= HYDRO_CUT else "polar")
            chg = ("exposed->buried" if exp_a[i] and not exp_b[i]
                   else "buried->exposed" if exp_b[i] and not exp_a[i]
                   else "exposed" if exp_a[i] and exp_b[i] else "buried")
            fh.write(f"{rn_a[i]},{names_a[i]},{charge[i]},{kd[i]:.1f},"
                     f"{sasa_a[i]:.2f},{sasa_b[i]:.2f},{sasa_b[i]-sasa_a[i]:+.2f},{cls},{chg}\n")

    # summary + patches
    out = [f"Surface charge / hydrophobic patch comparison: {tag}",
           f"  {labels[0]} = {structures[0]}",
           f"  {labels[1]} = {structures[1]}",
           f"  {n} residues; exposed if rSASA > {expose}; hydrophobic if KD >= {HYDRO_CUT}", ""]

    def section(title, ids_a, ids_b, coords_a, coords_b):
        out.append(f"== {title} ==")
        for lab, ids, co in ((labels[0], ids_a, coords_a), (labels[1], ids_b, coords_b)):
            res = ", ".join(rtag(i) for i in ids) or "(none)"
            out.append(f"  {lab}: {len(ids)} exposed [{res}]")
            for pi, patch in enumerate(find_patches(ids, co), 1):
                if len(patch) >= 2:
                    out.append(f"      patch {pi}: {'-'.join(rtag(i) for i in patch)}")
        gained = sorted(set(ids_b) - set(ids_a))
        lost = sorted(set(ids_a) - set(ids_b))
        out.append(f"  newly exposed in {labels[1]}: {', '.join(rtag(i) for i in gained) or '(none)'}")
        out.append(f"  buried on folding:          {', '.join(rtag(i) for i in lost) or '(none)'}")
        out.append("")

    section("Positive (Arg/Lys)", pos_a, pos_b, cb_a, cb_b)
    section("Negative (Asp/Glu)", neg_a, neg_b, cb_a, cb_b)
    section("Hydrophobic (KD>=1.5)", pho_a, pho_b, cb_a, cb_b)

    summ = os.path.join(outdir, f"{tag}_summary.txt")
    open(summ, "w").write("\n".join(out))
    if verbose:
        print("\n".join(out))
        print(f"\nwrote: {outdir}/{tag}_{{{','.join(labels)}}}.pdb, {tag}_residues.csv, {tag}_summary.txt")
    def patch_tags(ids, coords):          # patches = 2+ adjacent exposed residues of one type
        return [[rtag(i) for i in p] for p in find_patches(ids, coords) if len(p) >= 2]

    return {"n_res": n,
            "patches": {"apo": {"pos": patch_tags(pos_a, cb_a), "neg": patch_tags(neg_a, cb_a),
                                "hyd": patch_tags(pho_a, cb_a)},
                        "bound": {"pos": patch_tags(pos_b, cb_b), "neg": patch_tags(neg_b, cb_b),
                                  "hyd": patch_tags(pho_b, cb_b)}},
            "pos": [f"{names_a[i]}{rn_a[i]}" for i in pos_b],
            "neg": [f"{names_a[i]}{rn_a[i]}" for i in neg_b],
            "hydrophobic": [f"{names_a[i]}{rn_a[i]}" for i in pho_b],
            "buried_on_fold": [f"{names_a[i]}{rn_a[i]}" for i in range(n) if exp_a[i] and not exp_b[i]],
            "newly_exposed": [f"{names_a[i]}{rn_a[i]}" for i in range(n) if exp_b[i] and not exp_a[i]]}


def main():
    ap = argparse.ArgumentParser(description="Compare surface charge/hydrophobic patches between two peptide states.")
    ap.add_argument("structures", nargs=2, help="two states in order: apo cofold (.pdb or .cif)")
    ap.add_argument("--tag", default="peptide")
    ap.add_argument("--labels", default="apo,cofold", help="comma-separated state labels")
    ap.add_argument("--outdir", default="patches")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--expose", type=float, default=0.25, help="rSASA above which a residue is 'exposed'")
    args = ap.parse_args()
    analyze(args.structures, args.tag, args.labels.split(","), args.outdir, args.chain, args.expose)


if __name__ == "__main__":
    main()
