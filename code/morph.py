#!/usr/bin/env python3
"""
morph.py — build a PyMOL-viewable morph trajectory between peptide structures.

Takes an ordered list of structures (A B C ...), extracts the peptide chain
(heavy atoms only, ligand dropped), superposes each onto the first, and writes
one multi-MODEL PDB that linearly interpolates coordinates between consecutive
structures. Load the output in PyMOL and play it as a movie.

Handles both .pdb (apo folds) and Boltz-style .cif (co-folds). For .cif only the
polymer `ATOM` records are read, so the ligand `HETATM` block is dropped and you
get the peptide alone — exactly what we want here.

Atoms are matched across structures by (residue number, atom name) within the
chosen chain, so the morph is correct even if two programs wrote atoms in a
different order; only the atoms common to every structure are morphed.

Usage:
  morph.py -o traj.pdb --chain A --frames 30 apo.pdb cofold.cif
  morph.py -o traj.pdb --frames 24 esm.pdb of3.pdb rf3.pdb boltz.cif
Options:
  -o/--out     output multi-model PDB                 (required)
  --chain      chain to morph                          (default A)
  --frames     interpolation steps per segment         (default 30)
  --fit        superposition atoms: CA | backbone | all (default CA)
  --no-fit     skip superposition (use raw coordinates)
  --dwell      extra still frames held at each waypoint (default 0)
"""
import argparse
import sys
import numpy as np

BACKBONE = {"N", "CA", "C", "O"}


def parse_pdb(path, chain):
    """Return {(resseq:int, atomname:str): (name, resname, resseq, element, xyz)} for chain, ATOM only."""
    atoms = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):          # drops HETATM (ligands/ions) and headers
                continue
            if line[21] != chain:
                continue
            name = line[12:16].strip()
            resname = line[17:20].strip()
            resseq = int(line[22:26])
            element = line[76:78].strip() or name[0]
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            atoms[(resseq, name)] = (name, resname, resseq, element, xyz)
    return atoms


def parse_cif(path, chain):
    """Parse a Boltz/OpenFold-style mmCIF atom_site loop; peptide ATOM records only."""
    cols, rows, in_loop, reading = [], [], False, False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s == "loop_":
                in_loop, reading, cols = True, False, []
                continue
            if in_loop and s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1])
                continue
            if cols and (s.startswith("ATOM") or s.startswith("HETATM")):
                reading = True
                rows.append(s.split())
                continue
            if reading and (s == "#" or s.startswith("_") or s == "loop_" or not s):
                break  # end of the atom_site loop
    if not cols:
        sys.exit(f"{path}: no _atom_site loop found")
    idx = {c: i for i, c in enumerate(cols)}
    # auth_asym_id is the author chain id matching PDB; fall back to label_asym_id
    chain_col = idx.get("auth_asym_id", idx.get("label_asym_id"))
    seq_col = idx.get("auth_seq_id", idx.get("label_seq_id"))
    atoms = {}
    for r in rows:
        if r[idx["group_PDB"]] != "ATOM":            # drops HETATM ligand block
            continue
        if r[chain_col] != chain:
            continue
        name = r[idx["label_atom_id"]].strip('"')
        resname = r[idx["label_comp_id"]]
        resseq = int(r[seq_col])
        element = r[idx["type_symbol"]]
        xyz = (float(r[idx["Cartn_x"]]), float(r[idx["Cartn_y"]]), float(r[idx["Cartn_z"]]))
        atoms[(resseq, name)] = (name, resname, resseq, element, xyz)
    return atoms


def load(path, chain):
    return parse_cif(path, chain) if path.lower().endswith(".cif") else parse_pdb(path, chain)


def kabsch(P, Q):
    """Rotation R and translation t that best map P onto Q (Q ~ P@R.T + t)."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, _, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det((V @ Wt).T))
    D = np.diag([1.0, 1.0, d])
    R = (V @ D @ Wt).T
    return R, Q.mean(0) - R @ P.mean(0)


def pdb_line(serial, name, resname, resseq, element, xyz, chain):
    # strict PDB columns: name 13-16, resName 18-20, chain 22, resSeq 23-26,
    # x/y/z 31-54, occ 55-60, bfac 61-66, element 77-78
    nm = name if len(name) >= 4 else f" {name:<3}"   # 1-letter elements start in col 14
    return (f"ATOM  {serial:5d} {nm} {resname:>3} {chain}{resseq:4d}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2}\n")


def build_morph(structures, out, chain="A", frames=30, fit="CA", no_fit=False, dwell=0, verbose=True):
    """Write a multi-MODEL morph PDB between the given structures. Returns frame count."""
    if len(structures) < 2:
        raise ValueError("need at least two structures")
    structs = [load(p, chain) for p in structures]
    for p, s in zip(structures, structs):
        if not s:
            raise ValueError(f"{p}: no ATOM records on chain {chain}")

    common = set(structs[0])
    for s in structs[1:]:
        common &= set(s)
    keys = [k for k in structs[0] if k in common]   # keep first structure's native atom order
    if not keys:
        raise ValueError("no atoms shared across all structures")
    if verbose:
        for p, s in zip(structures, structs):
            extra = len(s) - len(keys)
            print(f"  {p}: {len(s)} atoms on chain {chain}"
                  + (f"  ({extra} not shared, dropped)" if extra else ""))
        print(f"morphing {len(keys)} atoms common to all {len(structs)} structures")

    coords = [np.array([s[k][4] for k in keys]) for s in structs]
    template = [structs[0][k] for k in keys]

    if not no_fit:
        if fit == "CA":
            mask = np.array([k[1] == "CA" for k in keys])
        elif fit == "backbone":
            mask = np.array([k[1] in BACKBONE for k in keys])
        else:
            mask = np.ones(len(keys), bool)
        for j in range(1, len(coords)):
            R, t = kabsch(coords[j][mask], coords[0][mask])
            coords[j] = coords[j] @ R.T + t

    out_frames = []
    for j in range(len(coords) - 1):
        A, B = coords[j], coords[j + 1]
        for _ in range(dwell):
            out_frames.append(A)
        start = 0 if j == 0 else 1                   # avoid duplicating the shared endpoint
        for i in range(start, frames + 1):
            out_frames.append(A + (i / frames) * (B - A))
    for _ in range(dwell):
        out_frames.append(coords[-1])

    with open(out, "w") as fh:
        fh.write("REMARK  morph trajectory generated by morph.py\n")
        fh.write(f"REMARK  waypoints: {' -> '.join(structures)}\n")
        for m, fr in enumerate(out_frames, 1):
            fh.write(f"MODEL     {m:4d}\n")
            for serial, (atom, xyz) in enumerate(zip(template, fr), 1):
                name, resname, resseq, element, _ = atom
                fh.write(pdb_line(serial, name, resname, resseq, element, xyz, chain))
            fh.write("ENDMDL\n")
    if verbose:
        print(f"wrote {len(out_frames)} frames to {out}")
    return len(out_frames)


def main():
    ap = argparse.ArgumentParser(description="Morph peptide structures into a PyMOL trajectory.")
    ap.add_argument("structures", nargs="+", help="structures in morph order (.pdb or .cif)")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--frames", type=int, default=30, help="steps per segment")
    ap.add_argument("--fit", choices=["CA", "backbone", "all"], default="CA")
    ap.add_argument("--no-fit", action="store_true")
    ap.add_argument("--dwell", type=int, default=0, help="still frames held at each waypoint")
    args = ap.parse_args()
    try:
        build_morph(args.structures, args.out, args.chain, args.frames, args.fit, args.no_fit, args.dwell)
    except ValueError as e:
        sys.exit(str(e))


if __name__ == "__main__":
    main()
