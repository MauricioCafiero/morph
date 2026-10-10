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
        present = [l for l in LANE_ORDER if l in lanes]
        lane_tuples = ", ".join(f'("{l}", "{LANE_COLORS[l]}", "{lanes[l]}")' for l in present)
        # Per-frame SS can't be refreshed live on a multi-state object (dss sets one ss per
        # atom and rebuild doesn't reliably redraw during playback). So we split each lane's
        # trajectory into one object per frame, bake dss into each (correct for that frame's
        # geometry), and animate by toggling visibility -- which plays reliably. The two real
        # endpoints (frames 1 and N, the apo and bound predictions) render exactly right; the
        # interpolated in-between frames are approximate (linear-morph geometry).
        out = [f"# {mol_pep}: apo lanes morphing into the Boltz co-fold ({n} states each).",
               f"# Run from this folder:  pymol {mol_pep}.pml   (opens paused; press play to morph)",
               f"# {len(present)} lanes: blue=ESMFold  salmon=OF3-apo  green=RF3-apo",
               "reinitialize",
               "python",
               "from pymol import cmd",
               "cmd.set('suspend_undo', 1)    # the movie toggles 123 objects per frame; recording that",
               "                              # as undo history fills memory and freezes playback/pause",
               f"n = {n}",
               f"lanes = [{lane_tuples}]        # (name, color, file); first is the alignment ref",
               "ref = lanes[0][0]",
               "for name, color, fn in lanes:",
               "    cmd.load(fn, '_' + name)",
               "for name, color, fn in lanes[1:]:          # overlay on the shared bound state",
               "    cmd.align('_' + name, '_' + ref, mobile_state=n, target_state=n)",
               "for name, color, fn in lanes:              # explode into per-frame objects",
               "    cmd.split_states('_' + name, prefix=name + '_f')",
               "    cmd.delete('_' + name)",
               "    for i in range(1, n + 1):",
               "        o = '%s_f%04d' % (name, i)",
               "        cmd.dss(o)                         # SS baked from THIS frame's geometry",
               "        cmd.hide('everything', o)",
               "        cmd.show('cartoon', o)",
               "        cmd.color(color, o)",
               "        cmd.disable(o)",
               "    cmd.group(name, name + '_f*')           # 41 frames tucked under ONE clickable name",
               "    cmd.disable(name)                       # lane off by default; click its name to show",
               "cmd.set('cartoon_fancy_helices', 1)",
               "cmd.bg_color('black')",
               "cmd.set('movie_fps', 10)",
               "cmd.mset('1 x%d' % n)",
               "for i in range(1, n + 1):",
               "    # show frame i of every lane; a lane whose group is off stays hidden, so clicking",
               "    # the group names is all you need -- the movie never touches the groups themselves",
               "    cmd.mdo(i, 'disable *_f*; enable *_f%04d' % i)",
               "cmd.enable('*_f0001')",
               "cmd.frame(1)",
               "cmd.orient()",
               "python end",
               "# All lanes start OFF. Click esmfold / of3apo / rf3apo in the object panel to show one,",
               "# then press play. Frame 1 = apo prediction, last frame = Boltz co-fold."]
        pml = os.path.join(MORPHS, f"{mol_pep}.pml")
        open(pml, "w").write("\n".join(out) + "\n")
        written.append((mol_pep, sorted(lanes), pml))

    print(f"wrote {len(written)} per-peptide PyMOL scripts to {MORPHS}")
    for mol_pep, lanes, _ in written:
        print(f"  {mol_pep}.pml  ({', '.join(lanes)})")


if __name__ == "__main__":
    main()
