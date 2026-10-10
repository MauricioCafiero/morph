#!/usr/bin/env python3
"""
run_batch.py — morph + strain + surface-patch analysis for every peptide that has
a Boltz co-fold, across all three apo lanes (ESMFold, OF3-apo, RF3-apo), side by side.

Discovers peptides as the intersection of a molecule's fold lanes and its Boltz
results, then for each (molecule, peptide, lane) writes under results/:

  results/<mol>/<peptide>/<lane>/
      morph.pdb            apo(lane) -> boltz co-fold trajectory (PyMOL movie)
      <peptide>_<lane>.pdb, <peptide>_boltz.pdb   patch-colored structures
      <peptide>_residues.csv, <peptide>_summary.txt
  results/strain_summary.csv     one row per (mol, peptide, lane): apo/cofold E, dE
  results/patch_summary.csv      exposed pos/neg/hydrophobic spots + buried-on-fold
  results/skipped.csv            peptides with no Boltz co-fold (e.g. BoltzGen arms)
  results/run.log                progress

Resumable: a (mol, peptide, lane) whose morph.pdb + summary already exist is skipped.
Run with the openmm-md env python (needs openmm, pdbfixer, mdtraj).
"""
import csv
import json
import os
import shutil
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import morph
import strain
import surface_patches as sp

# External data source (the peptidebuilder repo). Not part of this repo; override with
# MORPH_PEP_ROOT if peptidebuilder lives elsewhere.
PEP_ROOT = os.environ.get("MORPH_PEP_ROOT", os.path.expanduser("~/python_mac/peptidebuilder/runs"))
MOLS = ["octinoxate", "oxybenzone"]
LANES = ["esmfold", "of3apo", "rf3apo"]
RESULTS = os.path.join(REPO, "results")


def boltz_cif(mol, pep):
    p = f"{PEP_ROOT}/{mol}/boltz/boltz_results_{pep}/predictions/{pep}/{pep}_model_0.cif"
    return p if os.path.exists(p) else None


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(RESULTS, "run.log"), "a") as fh:
        fh.write(line + "\n")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    skipped = []
    done = 0

    jobs = []
    for mol in MOLS:
        fold_dir = f"{PEP_ROOT}/{mol}/folds"
        peps = sorted(p[:-4] for p in os.listdir(f"{fold_dir}/esmfold") if p.endswith(".pdb"))
        for pep in peps:
            cif = boltz_cif(mol, pep)
            if not cif:
                skipped.append((mol, pep, "no boltz co-fold (BoltzGen arm?)"))
                continue
            for lane in LANES:
                apo = f"{fold_dir}/{lane}/{pep}.pdb"
                if os.path.exists(apo):
                    jobs.append((mol, pep, lane, apo, cif))

    log(f"{len(jobs)} (mol,peptide,lane) jobs; {len(skipped)} peptides skipped (no co-fold)")
    morphs = os.path.join(RESULTS, "morphs")   # flat collection of trajectories for easy loading
    os.makedirs(morphs, exist_ok=True)
    boltz_e = {}   # (mol,pep) -> (E, rmsd): co-fold energy computed ONCE and reused across lanes

    def boltz_energy(mol, pep, cif):
        """Declash energy of the co-fold peptide, once per peptide (shared across its lanes)."""
        key = (mol, pep)
        if key not in boltz_e:
            bdir = os.path.join(RESULTS, mol, pep, "_boltz")
            os.makedirs(bdir, exist_ok=True)
            bpep = os.path.join(bdir, "boltz.pdb")
            strain.peptide_pdb(cif, "A", bpep)
            e, _, rmsd = strain.energy(bpep, 7.0, cache=os.path.join(bdir, "boltz_H.pdb"))
            boltz_e[key] = (e, rmsd)
        return boltz_e[key]

    for n, (mol, pep, lane, apo, cif) in enumerate(jobs, 1):
        outdir = os.path.join(RESULTS, mol, pep, lane)
        os.makedirs(outdir, exist_ok=True)
        metrics = os.path.join(outdir, "metrics.json")
        tag = f"{n}/{len(jobs)} {mol}/{pep}/{lane}"
        if os.path.exists(metrics):
            log(f"skip (done) {tag}")
            continue
        try:
            morph.build_morph([apo, cif], os.path.join(outdir, "morph.pdb"),
                              frames=30, fit="CA", dwell=5, verbose=False)
            shutil.copyfile(os.path.join(outdir, "morph.pdb"),
                            os.path.join(morphs, f"{mol}_{pep}_{lane}_to_boltz.pdb"))
            # strain: apo relaxed per lane; co-fold relaxed once per peptide and reused
            apo_pep = os.path.join(outdir, f"{lane}.pdb")
            strain.peptide_pdb(apo, "A", apo_pep)
            e_apo, _, rmsd_apo = strain.energy(apo_pep, 7.0, cache=os.path.join(outdir, f"{lane}_H.pdb"))
            e_boltz, rmsd_boltz = boltz_energy(mol, pep, cif)
            d_strain = e_boltz - e_apo
            pr = sp.analyze([apo, cif], tag=pep, labels=[lane, "boltz"], outdir=outdir, verbose=False)
            # the limited declash should barely move heavy atoms; flag any state that moved
            # past the guard (conformation may have shifted) or still ended up repulsive
            flags = [f"moved:{l}" for l, r in ((lane, rmsd_apo), ("boltz", rmsd_boltz)) if r > strain.MAX_HEAVY_RMSD]
            flags += [f"clash:{l}" for l, e in ((lane, e_apo), ("boltz", e_boltz)) if e > 0]
            rec = {"mol": mol, "peptide": pep, "lane": lane,
                   "E_apo_kcal": round(e_apo, 1), "E_boltz_kcal": round(e_boltz, 1),
                   "dE_strain_kcal": round(d_strain, 1),
                   "rmsd_apo": round(rmsd_apo, 2), "rmsd_boltz": round(rmsd_boltz, 2),
                   "note": ";".join(flags),
                   "n_pos": len(pr["pos"]), "n_neg": len(pr["neg"]), "n_hydrophobic": len(pr["hydrophobic"]),
                   "pos": pr["pos"], "neg": pr["neg"], "hydrophobic": pr["hydrophobic"],
                   "buried_on_fold": pr["buried_on_fold"], "newly_exposed": pr["newly_exposed"]}
            json.dump(rec, open(metrics, "w"), indent=0)   # written last = job's done-sentinel
            done += 1
            log(f"ok   {tag}  dE={d_strain:+.1f} kcal/mol  buried={len(pr['buried_on_fold'])}")
        except Exception as e:
            log(f"FAIL {tag}: {e!r}")
            traceback.print_exc()

    # rebuild aggregate CSVs from every metrics.json present (reap-proof: no dup rows)
    recs = []
    for mol in MOLS:
        base = os.path.join(RESULTS, mol)
        for pep in sorted(os.listdir(base)) if os.path.isdir(base) else []:
            for lane in LANES:
                m = os.path.join(base, pep, lane, "metrics.json")
                if os.path.exists(m):
                    recs.append(json.load(open(m)))
    with open(os.path.join(RESULTS, "strain_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mol", "peptide", "lane", "E_apo_kcal", "E_boltz_kcal", "dE_strain_kcal",
                    "rmsd_apo", "rmsd_boltz", "note"])
        for r in recs:
            w.writerow([r["mol"], r["peptide"], r["lane"], r["E_apo_kcal"], r["E_boltz_kcal"],
                        r["dE_strain_kcal"], r.get("rmsd_apo", ""), r.get("rmsd_boltz", ""), r.get("note", "")])
    with open(os.path.join(RESULTS, "patch_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mol", "peptide", "lane", "n_pos", "n_neg", "n_hydrophobic",
                    "buried_on_fold", "newly_exposed"])
        for r in recs:
            w.writerow([r["mol"], r["peptide"], r["lane"], r["n_pos"], r["n_neg"], r["n_hydrophobic"],
                        ";".join(r["buried_on_fold"]) or "-", ";".join(r["newly_exposed"]) or "-"])
    with open(os.path.join(RESULTS, "skipped.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mol", "peptide", "reason"])
        w.writerows(skipped)
    write_results_md(recs, skipped)
    log(f"DONE: {done} analyses this run; {len(recs)} total in results; skipped {len(skipped)} peptides (no co-fold)")


def load_kinds():
    """peptide -> design kind, from each molecule's fold manifest."""
    kinds = {}
    for mol in MOLS:
        man = f"{PEP_ROOT}/{mol}/folds/manifest.csv"
        if os.path.exists(man):
            for row in csv.DictReader(open(man)):
                kinds[(mol, row["peptide"])] = row.get("kind", "")
    return kinds


def write_results_md(recs, skipped):
    """One row per peptide: per-lane strain dE + a charge/hydrophobic patch description."""
    kinds = load_kinds()
    by_pep = {}
    for r in recs:
        by_pep.setdefault((r["mol"], r["peptide"]), {})[r["lane"]] = r

    lines = [
        "# Results — apo → Boltz co-fold strain and surface-patch changes", "",
        "One row per peptide. **Strain ΔE** = E(co-fold) − E(apo) in kcal/mol, ff14SB + GBn2,",
        "each predicted structure scored as-is with only a limited declash relaxation (heavy",
        "atoms move <0.6 Å; see README). Negative ⇒ the bound conformation is *lower* energy",
        "than the apo fold; positive ⇒ the apo fold must strain to reach the bound state.",
        "ΔE is given for each independent apo model (ESMFold / OF3-apo / RF3-apo).", "",
        "**Surface patches** are from static residue charge (Arg/Lys +, Asp/Glu −) and",
        "Kyte–Doolittle hydrophobicity gated by solvent exposure. 'Charged (co-fold)' lists the",
        "exposed charged residues in the bound state; 'Hydrophobic Δ on fold' is how many",
        "hydrophobic residues bury (−) or newly expose (+) going apo→co-fold, per lane.", "",
        "⚠ = a `note` flag in `strain_summary.csv` (clash or heavy-atom move past guard).", "",
        "| Mol | Peptide | Kind | ΔE esm | ΔE of3 | ΔE rf3 | Charged (co-fold) | Hydrophobic Δ on fold (esm/of3/rf3) |",
        "|-----|---------|------|-------:|-------:|-------:|-------------------|-------------------------------------|",
    ]

    def de_cell(r):
        if not r:
            return "—"
        v = r["dE_strain_kcal"]
        s = f"{v:+.0f}" if abs(v) >= 0.5 else "0"       # avoid "-0"
        return s + ("⚠" if r.get("note") else "")

    for (mol, pep) in sorted(by_pep):
        lanes = by_pep[(mol, pep)]
        any_r = next(iter(lanes.values()))
        pos, neg = any_r["pos"], any_r["neg"]          # co-fold exposed charges (lane-independent)
        charged = (f"{len(pos)}+ ({','.join(pos)})" if pos else "0+") + "; " + \
                  (f"{len(neg)}− ({','.join(neg)})" if neg else "0−")
        hyd = []
        for lane in LANES:
            r = lanes.get(lane)
            if not r:
                hyd.append("—")
            else:
                nb, ne = len(r["buried_on_fold"]), len(r["newly_exposed"])
                hyd.append(f"−{nb}/+{ne}")
        lines.append(f"| {mol} | {pep} | {kinds.get((mol, pep), '')} | "
                     f"{de_cell(lanes.get('esmfold'))} | {de_cell(lanes.get('of3apo'))} | "
                     f"{de_cell(lanes.get('rf3apo'))} | {charged} | {' / '.join(hyd)} |")

    if skipped:
        lines += ["", "## Skipped (no Boltz co-fold — BoltzGen arms)", "",
                  ", ".join(f"{m}/{p}" for m, p, _ in skipped)]
    lines += ["", "---", "",
              "Per-peptide detail (patch residue lists, per-lane CSV, colored PDBs, morph movie)",
              "lives in `results/<mol>/<peptide>/<lane>/`; all morph trajectories are collected",
              "in `results/morphs/`. Full tables: `strain_summary.csv`, `patch_summary.csv`."]
    open(os.path.join(RESULTS, "RESULTS.md"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
