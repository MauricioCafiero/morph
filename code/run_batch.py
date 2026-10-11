#!/usr/bin/env python3
"""
run_batch.py — morph + strain + surface-patch analysis for every peptide, pairing each
apo fold with its own model's co-fold (Boltz apo -> Boltz, OF3 apo -> OF3, RF3 apo -> RF3).

The bound reference is that pair's co-fold; for BoltzGen peptides (bg33_N, bgox31_N) in
the Boltz pair it is BoltzGen's own co-fold from boltzgen_local, no Boltz co-fold existing
for those sequences. `--set legacy` reruns the superseded first analysis, which morphed
ESMFold / OF3-apo / RF3-apo all to the Boltz co-fold.
Input roots can be overridden with MORPH_PEP_ROOT and MORPH_BOLTZGEN_ROOT. For each
(molecule, peptide, lane) it writes under results/:

  results/<mol>/<peptide>/<lane>/
      morph.pdb            apo(lane) -> bound trajectory (PyMOL movie)
      <peptide>_<lane>.pdb, <peptide>_boltz.pdb   patch-colored structures
      <peptide>_residues.csv, <peptide>_summary.txt
      <lane>.pdb, <lane>_H.pdb                     peptide + cached protonation (strain)
      metrics.json                                 all numbers; written last = done-sentinel
  results/<mol>/<peptide>/_boltz/  bound structure relaxed once, shared by all lanes
  results/morphs/                  every morph movie + one PyMOL script per peptide
  results/strain_summary.csv       per (mol, peptide, lane): E_apo, E_bound, dE, RMSDs, flags
  results/patch_summary.csv        exposed +/-/hydrophobic counts, buried / newly exposed
  results/RESULTS.md               one row per peptide: strain per lane + Apo/Bound patches
  results/skipped.csv              peptides with no bound structure found
  results/run.log                  progress

Resumable: jobs with a metrics.json are skipped, and the CSVs / RESULTS.md / PyMOL
scripts are rebuilt from all metrics.json present. Run with the openmm-md env python
(needs openmm, pdbfixer, mdtraj).
"""
import argparse
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
LANES = ["esmfold", "of3apo", "rf3apo"]          # the apo lanes of the legacy set

# A *set* is a list of (apo lane, bound lane) pairs. The original run had no Boltz apo fold,
# so all three apo lanes were morphed to the one Boltz co-fold -- which compares an apo fold
# from one model against a bound fold from another. Since 2026-10-11 every co-folder has its
# own apo lane, so the matched set pairs each apo fold with ITS OWN model's co-fold and the
# apo->bound difference is the ligand's effect rather than a change of model. ESMFold has no
# co-fold, so it appears only in the legacy set.
LEGACY_PAIRS = [("esmfold", "boltz"), ("of3apo", "boltz"), ("rf3apo", "boltz")]
MATCHED_PAIRS = [("boltzapo", "boltz"), ("of3apo", "openfold3"), ("rf3apo", "rf3")]
SETS = {"legacy": LEGACY_PAIRS, "matched": MATCHED_PAIRS}

# both set by main() from --set; module level so the helpers below can see them
PAIRS = MATCHED_PAIRS
SETNAME = "matched"
RESULTS = os.path.join(REPO, "results")


# BoltzGen peptides have no Boltz co-fold; their bound reference is BoltzGen's own design
# co-fold (bg33_N -> octinoxate33_N.cif, bgox31_N -> oxybenzone31_N.cif), matched to the apo
# folds by sequence. Override with MORPH_BOLTZGEN_ROOT if boltzgen_local lives elsewhere.
BOLTZGEN_ROOT = os.environ.get("MORPH_BOLTZGEN_ROOT",
                               os.path.expanduser("~/python_mac/boltzgen_local/results"))
BOLTZGEN_PREFIX = {"octinoxate": ("bg33_", "octinoxate33_"), "oxybenzone": ("bgox31_", "oxybenzone31_")}


def boltz_cif(mol, pep):
    """Bound reference: the Boltz co-fold, or BoltzGen's co-fold for BoltzGen peptides."""
    p = f"{PEP_ROOT}/{mol}/boltz/boltz_results_{pep}/predictions/{pep}/{pep}_model_0.cif"
    if os.path.exists(p):
        return p
    pre, cif_pre = BOLTZGEN_PREFIX.get(mol, (None, None))
    if pre and pep.startswith(pre):
        p = f"{BOLTZGEN_ROOT}/{cif_pre}{pep[len(pre):]}.cif"
        if os.path.exists(p):
            return p
    return None


def bound_struct(mol, pep, bound):
    """Bound reference for one (peptide, bound lane): the Boltz co-fold (with BoltzGen's own
    co-fold standing in for BoltzGen peptides), or a cofolder lane's converted structure."""
    if bound == "boltz":
        return boltz_cif(mol, pep)
    p = f"{PEP_ROOT}/{mol}/folds/{bound}/{pep}.pdb"
    return p if os.path.exists(p) else None


def job_dirname(apo, bound):
    """Legacy job dirs are named by the apo lane alone; matched dirs name both ends, so the
    two sets can never collide even if they are written into the same tree."""
    return apo if SETNAME == "legacy" else f"{apo}_to_{bound}"


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(RESULTS, "run.log"), "a") as fh:
        fh.write(line + "\n")


def main():
    """Find peptides, run every missing (mol, peptide, pair) job, then rebuild all summaries."""
    global PAIRS, SETNAME, RESULTS
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", choices=sorted(SETS), default="matched", dest="setname",
                    help="matched (default): each apo lane -> its own model's co-fold, written "
                         "to results/. legacy: all apo lanes -> the Boltz co-fold, written to "
                         "results/legacy/ -- the superseded original run, kept runnable")
    ap.add_argument("--results", help="override the output directory")
    args = ap.parse_args()
    SETNAME = args.setname
    PAIRS = SETS[SETNAME]
    RESULTS = args.results or (os.path.join(REPO, "results") if SETNAME == "matched"
                               else os.path.join(REPO, "results", SETNAME))
    os.makedirs(RESULTS, exist_ok=True)
    skipped = []
    done = 0

    jobs = []
    for mol in MOLS:
        fold_dir = f"{PEP_ROOT}/{mol}/folds"
        peps = sorted(p[:-4] for p in os.listdir(f"{fold_dir}/{PAIRS[0][0]}") if p.endswith(".pdb"))
        for pep in peps:
            for lane, bound in PAIRS:
                apo = f"{fold_dir}/{lane}/{pep}.pdb"
                ref = bound_struct(mol, pep, bound)
                if not os.path.exists(apo):
                    continue
                if not ref:
                    skipped.append((mol, pep, f"no {bound} bound structure for the {lane} pair"))
                    continue
                jobs.append((mol, pep, lane, bound, apo, ref))

    log(f"set={SETNAME}: {len(jobs)} (mol,peptide,pair) jobs; {len(skipped)} pairs skipped")
    morphs = os.path.join(RESULTS, "morphs")   # flat collection of trajectories for easy loading
    os.makedirs(morphs, exist_ok=True)
    boltz_e = {}   # (mol,pep,bound) -> (E, rmsd): each bound structure relaxed ONCE

    def boltz_energy(mol, pep, bound, ref):
        """Declash energy of a bound peptide, once per (peptide, bound lane)."""
        key = (mol, pep, bound)
        if key not in boltz_e:
            bdir = os.path.join(RESULTS, mol, pep, f"_{bound}")
            os.makedirs(bdir, exist_ok=True)
            bpep = os.path.join(bdir, f"{bound}.pdb")
            strain.peptide_pdb(ref, "A", bpep)
            e, _, rmsd = strain.energy(bpep, 7.0, cache=os.path.join(bdir, f"{bound}_H.pdb"))
            boltz_e[key] = (e, rmsd)
        return boltz_e[key]

    for n, (mol, pep, lane, bound, apo, cif) in enumerate(jobs, 1):
        outdir = os.path.join(RESULTS, mol, pep, job_dirname(lane, bound))
        os.makedirs(outdir, exist_ok=True)
        metrics = os.path.join(outdir, "metrics.json")
        tag = f"{n}/{len(jobs)} {mol}/{pep}/{lane}->{bound}"
        if os.path.exists(metrics):
            log(f"skip (done) {tag}")
            continue
        try:
            morph.build_morph([apo, cif], os.path.join(outdir, "morph.pdb"),
                              frames=30, fit="CA", dwell=5, verbose=False)
            shutil.copyfile(os.path.join(outdir, "morph.pdb"),
                            os.path.join(morphs, f"{mol}_{pep}_{lane}_to_{bound}.pdb"))
            # strain: apo relaxed per lane; co-fold relaxed once per peptide and reused
            apo_pep = os.path.join(outdir, f"{lane}.pdb")
            strain.peptide_pdb(apo, "A", apo_pep)
            e_apo, _, rmsd_apo = strain.energy(apo_pep, 7.0, cache=os.path.join(outdir, f"{lane}_H.pdb"))
            e_boltz, rmsd_boltz = boltz_energy(mol, pep, bound, cif)
            d_strain = e_boltz - e_apo
            pr = sp.analyze([apo, cif], tag=pep, labels=[lane, bound], outdir=outdir, verbose=False)
            # the limited declash should barely move heavy atoms; flag any state that moved
            # past the guard (conformation may have shifted) or still ended up repulsive
            flags = [f"moved:{l}" for l, r in ((lane, rmsd_apo), (bound, rmsd_boltz)) if r > strain.MAX_HEAVY_RMSD]
            flags += [f"clash:{l}" for l, e in ((lane, e_apo), (bound, e_boltz)) if e > 0]
            rec = {"mol": mol, "peptide": pep, "lane": lane, "bound": bound,
                   "E_apo_kcal": round(e_apo, 1), "E_boltz_kcal": round(e_boltz, 1),
                   "dE_strain_kcal": round(d_strain, 1),
                   "rmsd_apo": round(rmsd_apo, 2), "rmsd_boltz": round(rmsd_boltz, 2),
                   "note": ";".join(flags),
                   "n_pos": len(pr["pos"]), "n_neg": len(pr["neg"]), "n_hydrophobic": len(pr["hydrophobic"]),
                   "patches": pr["patches"],
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
            for lane, bound in PAIRS:
                m = os.path.join(base, pep, job_dirname(lane, bound), "metrics.json")
                if os.path.exists(m):
                    recs.append(json.load(open(m)))
    with open(os.path.join(RESULTS, "strain_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mol", "peptide", "lane", "bound", "E_apo_kcal", "E_bound_kcal", "dE_strain_kcal",
                    "rmsd_apo", "rmsd_bound", "note"])
        for r in recs:
            w.writerow([r["mol"], r["peptide"], r["lane"], r.get("bound", "boltz"), r["E_apo_kcal"],
                        r["E_boltz_kcal"], r["dE_strain_kcal"], r.get("rmsd_apo", ""),
                        r.get("rmsd_boltz", ""), r.get("note", "")])
    with open(os.path.join(RESULTS, "patch_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mol", "peptide", "lane", "bound", "n_pos", "n_neg", "n_hydrophobic",
                    "buried_on_fold", "newly_exposed"])
        for r in recs:
            w.writerow([r["mol"], r["peptide"], r["lane"], r.get("bound", "boltz"),
                        r["n_pos"], r["n_neg"], r["n_hydrophobic"],
                        ";".join(r["buried_on_fold"]) or "-", ";".join(r["newly_exposed"]) or "-"])
    with open(os.path.join(RESULTS, "skipped.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mol", "peptide", "reason"])
        w.writerows(skipped)
    write_results_md(recs, skipped)
    try:
        import make_pymol
        # this set's morphs dir only -- never another set's scripts
        make_pymol.main(morphs=morphs, pairs=PAIRS)
    except SystemExit:
        pass
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

    SHORT = {"esmfold": "esm", "of3apo": "of3", "rf3apo": "rf3", "boltzapo": "boltz",
             "boltz": "boltz", "openfold3": "of3", "rf3": "rf3"}
    # matched pairs live inside one model, so the column names it once; the legacy set's
    # columns name the apo lane, the bound end being the same Boltz co-fold throughout
    cols = " | ".join(f"ΔE {SHORT[a]}" if SETNAME == "matched" else f"ΔE {SHORT[a]}"
                      for a, b in PAIRS)
    rule = "|".join("-------:" for _ in PAIRS)
    if SETNAME == "matched":
        intro = [
            "# Results — matched apo → co-fold strain and surface-patch changes", "",
            "Each column pairs an apo fold with **its own model's** co-fold: Boltz apo → Boltz",
            "co-fold, OF3 apo → OF3 co-fold, RF3 apo → RF3 co-fold. The apo→bound difference is",
            "therefore the ligand's effect within one model, not a change of model. The legacy",
            "set in `results/` morphed all three apo lanes to the single Boltz co-fold, which",
            "mixed the two; it is kept for comparison. ESMFold has no co-fold, so it is absent",
            "here and present only in the legacy set.", "",
            "For the Boltz pair, BoltzGen peptides (kind = boltzgen: bg33_*, bgox31_*) keep",
            "BoltzGen's own design co-fold as the bound reference — no Boltz co-fold exists for",
            "those sequences — so that one column is matched for 13 of 21 peptides and",
            "BoltzGen-bound for the other 8.", ""]
    else:
        intro = ["# Results — apo → Boltz co-fold strain and surface-patch changes", ""]
    lines = intro + [
        "One row per peptide. **Strain ΔE** = E(co-fold) − E(apo) in kcal/mol, ff14SB + GBn2,",
        "each predicted structure scored as-is with only a limited declash relaxation (heavy",
        "atoms move <0.6 Å; see README). Negative ⇒ the bound conformation is *lower* energy",
        "than the apo fold; positive ⇒ the apo fold must strain to reach the bound state.",
        "ΔE is given for each (apo model → bound model) pair in this set.", "",
        "**Bound reference.** The Boltz co-fold for designed / variant / null peptides; for",
        "BoltzGen peptides (kind = boltzgen: bg33_*, bgox31_*) it is BoltzGen's own design",
        "co-fold from `boltzgen_local/results/` (matched to the apo folds by sequence). Columns",
        "and files still say \"boltz\" for both.", "",
        "**Surface patches.** A patch is 2+ solvent-exposed residues of one type whose Cβ atoms",
        "lie within 8 Å of each other: **+** (Arg/Lys), **−** (Asp/Glu), or **hydrophobic**",
        "(Kyte–Doolittle ≥ 1.5). Isolated exposed residues are not counted. **Apo** gives the",
        "patch counts on each apo model's surface. **Bound** says what happened to them in the",
        "Boltz co-fold: *buried* (most of the patch's residues are no longer exposed), *still",
        "exposed* (the patch persists), *broken up* (residues still exposed but no longer",
        "clustered), and *new* (a patch that only exists in the bound state). The bound",
        ("structure is the same for all three lanes; the changes differ because the apo models do."
         if SETNAME == "legacy" else
         "bound structure differs per pair, so both ends of each comparison come from one model."), "",
        "⚠ = a `note` flag in `strain_summary.csv` (clash or heavy-atom move past guard).", "",
        f"| Mol | Peptide | Kind | {cols} | Apo (patches) | Bound (change on folding) |",
        f"|-----|---------|------|{rule}|---------------|---------------------------|",
    ]
    LABEL = {"pos": "+", "neg": "−", "hyd": "hydrophobic"}

    def plural(k, word="patch"):
        return f"{k} {word}" + ("es" if k != 1 else "")

    def apo_cell(r):
        a = r["patches"]["apo"]
        return f"{len(a['pos'])} +, {len(a['neg'])} −, {len(a['hyd'])} hydrophobic"

    def bound_cell(r):
        """Fate of each apo patch in the bound state: >=50% of its residues in one bound patch
        -> still exposed; else >=50% buried -> buried; else broken up. Bound patches with
        <50% overlap with every apo patch -> new."""
        a, b = r["patches"]["apo"], r["patches"]["bound"]
        buried = set(r["buried_on_fold"])
        parts = []
        for c in ("pos", "neg", "hyd"):
            fate = {"buried": 0, "still exposed": 0, "broken up": 0, "new": 0}
            for p in a[c]:
                keep = max((len(set(p) & set(q)) / len(p) for q in b[c]), default=0)
                if keep >= 0.5:
                    fate["still exposed"] += 1
                elif len(set(p) & buried) / len(p) >= 0.5:
                    fate["buried"] += 1
                else:
                    fate["broken up"] += 1
            for q in b[c]:
                if max((len(set(q) & set(p)) / len(q) for p in a[c]), default=0) < 0.5:
                    fate["new"] += 1
            bits = [f"{k} {v}" for v, k in fate.items() if k]
            if bits:
                parts.append(f"{LABEL[c]}: " + ", ".join(bits))
        return "; ".join(parts) or "no patches"

    def per_lane(lanes, fn):
        return "<br>".join(f"**{SHORT[a]}** {fn(lanes[a])}" for a, _ in PAIRS if a in lanes)

    def de_cell(r):
        if not r:
            return "—"
        v = r["dE_strain_kcal"]
        s = f"{v:+.0f}" if abs(v) >= 0.5 else "0"       # avoid "-0"
        return s + ("⚠" if r.get("note") else "")

    for (mol, pep) in sorted(by_pep):
        lanes = by_pep[(mol, pep)]
        cells = " | ".join(de_cell(lanes.get(a)) for a, _ in PAIRS)
        lines.append(f"| {mol} | {pep} | {kinds.get((mol, pep), '')} | {cells} | "
                     f"{per_lane(lanes, apo_cell)} | {per_lane(lanes, bound_cell)} |")

    if skipped:
        lines += ["", "## Skipped (no bound structure found)", "",
                  ", ".join(f"{m}/{p}" for m, p, _ in skipped)]
    lines += ["", "---", "",
              "Per-peptide detail (patch residue lists, per-lane CSV, colored PDBs, morph movie)",
              "lives in `results/<mol>/<peptide>/<lane>/`; all morph trajectories are collected",
              "in `results/morphs/`. Full tables: `strain_summary.csv`, `patch_summary.csv`."]
    open(os.path.join(RESULTS, "RESULTS.md"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
