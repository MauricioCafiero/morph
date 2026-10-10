#!/usr/bin/env python3
"""
make_pymol.py — write PyMOL scripts that load each peptide's apo-lane morphs together.

Scans results/morphs/ (files named <mol>_<peptide>_<lane>_to_boltz.pdb) and, for every
peptide, writes results/morphs/<mol>_<peptide>.pml which loads that peptide's three lane
morphs (ESMFold / OF3-apo / RF3-apo) as separate multi-state objects, brings them into a
common frame via their shared final (Boltz-bound) state, colors them per lane, and leaves
them ready to play as a synchronized movie.

Run from the repo root (reads what run_batch.py produced):
    python code/make_pymol.py
Then, to view one peptide (from the morphs folder so the bare filenames resolve):
    cd results/morphs && pymol octinoxate_s3_esm2_f4.pml      # press play
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MORPHS = os.path.join(os.path.dirname(HERE), "results", "morphs")

LANE_COLORS = {"esmfold": "skyblue", "of3apo": "salmon", "rf3apo": "palegreen"}
LANE_ORDER = ["esmfold", "of3apo", "rf3apo"]


def n_states(pdb_path):
    with open(pdb_path) as fh:
        return sum(1 for line in fh if line.startswith("MODEL"))


def main():
    if not os.path.isdir(MORPHS):
        sys.exit(f"no morphs folder at {MORPHS} — run code/run_batch.py first")

    # group morph files by (mol, peptide) -> {lane: filename}
    peptides = {}
    for fn in sorted(os.listdir(MORPHS)):
        if not fn.endswith("_to_boltz.pdb"):
            continue
        stem = fn[:-len("_to_boltz.pdb")]
        for lane in LANE_ORDER:
            if stem.endswith("_" + lane):
                mol_pep = stem[:-(len(lane) + 1)]
                peptides.setdefault(mol_pep, {})[lane] = fn
                break

    if not peptides:
        sys.exit(f"no '*_to_boltz.pdb' morphs found in {MORPHS}")

    written = []
    for mol_pep, lanes in sorted(peptides.items()):
        ref_lane = next(l for l in LANE_ORDER if l in lanes)
        n = n_states(os.path.join(MORPHS, lanes[ref_lane]))
        out = [f"# {mol_pep}: apo lanes morphing into the Boltz co-fold ({n} states each).",
               f"# Run from this folder:  pymol {mol_pep}.pml   then press the play button.",
               "reinitialize", ""]
        for lane in LANE_ORDER:
            if lane in lanes:
                out.append(f"load {lanes[lane]}, {lane}")
        out.append("")
        # overlay all lanes on their shared final (bound) state, so they share one frame
        for lane in LANE_ORDER:
            if lane in lanes and lane != ref_lane:
                out.append(f"align {lane}, {ref_lane}, mobile_state={n}, target_state={n}")
        out += ["", "hide everything", "show cartoon",
                "dss                          # assign secondary structure (PDBs carry none) so",
                "                             # helices/strands render as ribbon, not thin coil",
                "set cartoon_transparency, 0",
                "set cartoon_fancy_helices, 1"]
        for lane in LANE_ORDER:
            if lane in lanes:
                out.append(f"color {LANE_COLORS[lane]}, {lane}")
        out += ["", "set all_states, off          # show one movie frame at a time",
                "bg_color black", "set ray_opaque_background, 0", "set movie_fps, 15", "orient",
                "set scene_buttons, 1",
                "# the three lanes share states 1..%d; drag the slider or press play to morph" % n]
        pml = os.path.join(MORPHS, f"{mol_pep}.pml")
        open(pml, "w").write("\n".join(out) + "\n")
        written.append((mol_pep, sorted(lanes), pml))

    print(f"wrote {len(written)} per-peptide PyMOL scripts to {MORPHS}")
    for mol_pep, lanes, _ in written:
        print(f"  {mol_pep}.pml  ({', '.join(lanes)})")


if __name__ == "__main__":
    main()
