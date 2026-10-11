#!/usr/bin/env python3
"""
make_pymol.py — write one PyMOL script per peptide for its three apo -> bound morphs,
plus one patch viewer per lane job directory.

Scans results/morphs/ (files named <mol>_<peptide>_<lane>_to_boltz.pdb) and, for every
peptide, writes results/morphs/<mol>_<peptide>.pml. The script:
  * loads the three lane morphs (ESMFold / OF3-apo / RF3-apo) and aligns them on their
    shared final (bound) state;
  * splits each lane into one object per frame and runs dss on each, so secondary
    structure is computed per frame (PyMOL stores one ss per atom, so a single
    multi-state object can't carry per-frame secondary structure);
  * groups each lane's frame objects under one clickable name (esmfold / of3apo / rf3apo),
    all off by default;
  * builds a movie that shows frame i of each lane and hides its neighbours, so a lane
    whose group is off stays hidden and clicking group names is the whole interface;
  * opens paused on frame 1. Play with mplay; stop with mstop (the GUI pause button
    doesn't stop it). suspend_undo is set where the PyMOL build supports it.

Run from the repo root (reads what run_batch.py produced):
    python code/make_pymol.py
Then view one peptide from the morphs folder so the bare filenames resolve:
    cd results/morphs && pymol octinoxate_s3_esm2_f4.pml

It also writes <peptide>_patches.pml into each results/<mol>/<peptide>/<lane>/ job
directory (the loader for that pair of patch-colored PDBs, see surface_patches.py) so
already-run jobs always have it. Then view the patch map from that job's folder:
    cd results/oxybenzone/ox2_esm1_f8/of3apo && pymol ox2_esm1_f8_patches.pml
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# The morphs directory of ONE result set. run_batch passes its own via main(morphs=...);
# MORPH_MORPHS_DIR overrides it for a standalone run. Defaults to the legacy set.
MORPHS = os.environ.get("MORPH_MORPHS_DIR",
                        os.path.join(os.path.dirname(HERE), "results", "morphs"))

LANE_COLORS = {"esmfold": "skyblue", "of3apo": "salmon", "rf3apo": "palegreen",
               "boltzapo": "wheat"}
LANE_ORDER = ["esmfold", "of3apo", "rf3apo", "boltzapo"]
# bound lanes a morph file can end in: <mol>_<pep>_<apo lane>_to_<bound lane>.pdb
BOUND_LANES = ["boltz", "openfold3", "rf3"]


def n_states(pdb_path):
    with open(pdb_path) as fh:
        return sum(1 for line in fh if line.startswith("MODEL"))


def main(morphs=None, pairs=None):
    """Write <morphs>/<mol>_<peptide>.pml for every peptide with morphs there.

    `morphs` is one set's morphs directory (default: the module-level MORPHS, i.e. the legacy
    set) and `pairs` the (apo, bound) pairs that set holds, used to find each job directory
    for the patch viewers. Writing is confined to the directory passed, so running one set
    never rewrites another's scripts."""
    MORPHS = morphs or globals()["MORPHS"]
    pairs = pairs or [(l, "boltz") for l in ("esmfold", "of3apo", "rf3apo")]
    if not os.path.isdir(MORPHS):
        sys.exit(f"no morphs folder at {MORPHS} — run code/run_batch.py first")

    # group morph files by (mol, peptide) -> {apo lane: filename}, over any bound lane
    peptides = {}
    for fn in sorted(os.listdir(MORPHS)):
        hit = next((b for b in BOUND_LANES if fn.endswith(f"_to_{b}.pdb")), None)
        if not hit:
            continue
        stem = fn[:-len(f"_to_{hit}.pdb")]
        for lane in LANE_ORDER:
            if stem.endswith("_" + lane):
                mol_pep = stem[:-(len(lane) + 1)]
                peptides.setdefault(mol_pep, {})[lane] = fn
                break

    if not peptides:
        sys.exit(f"no '*_to_<bound>.pdb' morphs found in {MORPHS}")

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
               "try:                          # undo history of per-frame toggles fills memory and",
               "    cmd.set('suspend_undo', 1) # stalls pause; not every PyMOL build has this setting",
               "except Exception:",
               "    pass",
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
               "cmd.set('movie_panel', 1)                 # timeline with play/stop controls",
               "cmd.mset('1 x%d' % n)",
               "for i in range(1, n + 1):",
               "    # show frame i of every lane and hide only its neighbours (few toggles per frame,",
               "    # so undo history grows slowly); frame 1 also clears everything once per loop.",
               "    # A lane whose group is off stays hidden, so clicking group names is the interface.",
               "    prev = n if i == 1 else i - 1",
               "    nxt = 1 if i == n else i + 1",
               "    reset = 'disable *_f*; ' if i == 1 else ''",
               "    cmd.mdo(i, reset + 'disable *_f%04d; disable *_f%04d; enable *_f%04d' % (prev, nxt, i))",
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

    # patch viewers: one per lane job dir that has the patch-colored pair (surface_patches
    # writes them when a job runs; this pass covers jobs done earlier)
    import surface_patches as sp
    RESULTS = os.path.dirname(MORPHS)
    n_patches = 0
    for mol in sorted(os.listdir(RESULTS)):
        mol_dir = os.path.join(RESULTS, mol)
        if not os.path.isdir(mol_dir):
            continue
        for pep in sorted(os.listdir(mol_dir)):
            if pep.startswith("_"):
                continue  # shared _boltz dir
            pep_dir = os.path.join(mol_dir, pep)
            if not os.path.isdir(pep_dir):
                continue
            for lane, bound in pairs:
                # legacy job dirs are named by the apo lane; matched ones by both ends
                for dirname in (lane, f"{lane}_to_{bound}"):
                    lane_dir = os.path.join(pep_dir, dirname)
                    if os.path.isfile(os.path.join(lane_dir, f"{pep}_{lane}.pdb")):
                        sp.write_patch_pml(pep, [lane, bound], lane_dir)
                        n_patches += 1
    if n_patches:
        print(f"wrote {n_patches} patch viewer scripts into results/<mol>/<pep>/<lane>/")
    for mol_pep, lanes, _ in written:
        print(f"  {mol_pep}.pml  ({', '.join(lanes)})")


if __name__ == "__main__":
    main()
