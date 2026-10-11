# morph — apo → co-fold peptide analysis

For each designed peptide, compare an **apo** prediction with the peptide's **bound** structure
from a co-fold with its ligand, and produce:

1. **Morph movies** — apo → bound trajectories you can play in PyMOL, with secondary structure
   computed per frame.
2. **Strain** — the MM energy cost (or gain) for the peptide to go from its apo shape to the
   bound shape (ff14SB + GBn2).
3. **Surface patches** — charged (+/−) and hydrophobic patches on the apo surface, and whether
   they are buried, persist, break up, or newly form in the bound state.

There are two sets of pairs, selected with `--set`:

| set | pairs | output | what a ΔE means |
|---|---|---|---|
| `matched` (default) | Boltz-apo → Boltz, OF3-apo → OF3, RF3-apo → RF3 | `results/` | the ligand's effect **within one model** |
| `legacy` | ESMFold, OF3-apo, RF3-apo → **the Boltz co-fold** | `results/legacy/` | superseded: mixes the ligand's effect with a change of model between the apo and bound ends |

The legacy set came first, when no Boltz apo fold existed: all three apo lanes were morphed to
the single Boltz co-fold, so an apo fold from one model was compared against a bound fold from
another. Every co-folder now has its own apo lane (Boltz apo added 2026-10-11), so the matched
set pairs each apo fold with its own model's co-fold and the difference is attributable to the
ligand. **The matched set replaces it** (2026-10-11).

The legacy set's distinguishing feature was being the only place ESMFold appears, having no
co-fold of its own — and that is now a reason to retire it rather than keep it. Inverse-folding
perplexity over all six fold lanes (the `pepgraph` repo) puts ESMFold last by a wide margin: it
owns 7 of the 11 structures whose geometry carries no sequence information at all, and 6 of
those have *exactly zero* residue contacts, Rg 19–28 Å. It emits extended chains without
secondary structure for these ligand-dependent designs, so morphing from one is measuring a
model failure, not an apo state. The set stays runnable for provenance, and its outputs are
regenerable with `--set legacy`, but nothing should be concluded from it.

In the matched set's Boltz column, BoltzGen peptides (`bg33_*`, `bgox31_*`) keep BoltzGen's own
design co-fold as the bound reference, there being no Boltz co-fold for those sequences — so
that column is model-matched for 13 of 21 peptides and BoltzGen-bound for the other 8.

Each set is summarized one row per peptide in its own `RESULTS.md` (`results/RESULTS.md` for
the matched set, `results/legacy/RESULTS.md` for the superseded one). Strain is meant as a
qualitative cross-check on the binding energetics computed elsewhere (MM/GBSA in
`peptidebuilder`, absolute binding free energy in `fep`), not as a binding energy itself.

---

## Workflow

### 1. Inputs (read in place, never modified)

The inputs are prediction outputs from separate local project trees (peptide fold runs,
Boltz / BoltzGen co-folds); they are not shipped in this repo. Point the scripts at them
with two environment variables:

| What | Layout relative to its root | Root variable |
|---|---|---|
| Apo folds, one PDB per peptide per lane | `<mol>/folds/{esmfold,of3apo,rf3apo,boltzapo}/<peptide>.pdb` | `MORPH_PEP_ROOT` (default `<home>/<project>/runs`) |
| Co-fold lanes for the matched set | `<mol>/folds/{openfold3,rf3}/<peptide>.pdb` (ligand on chain Z / B, dropped) | same |
| Boltz co-folds (designed / variant / null peptides) | `<mol>/boltz/boltz_results_<peptide>/predictions/<peptide>/<peptide>_model_0.cif` | same |
| BoltzGen co-folds (`bg33_N`, `bgox31_N`) | `<name>_N.cif` in a flat result directory | `MORPH_BOLTZGEN_ROOT` |
| Peptide "kind" for the results table | `<mol>/folds/manifest.csv` | — |

`<mol>` is `octinoxate` or `oxybenzone`. Only the peptide chain (A) is used; the ligand is
dropped. BoltzGen peptides have no Boltz co-fold, so their bound reference is BoltzGen's own
design co-fold (matched to the apo folds by sequence). Files and columns still say "boltz"
for both.

Concrete paths on the current machine are in `LOCAL.md` (untracked).

### 2. Environment

Needs a python environment with `openmm`, `pdbfixer`, `mdtraj` and `numpy`. Set `$PY` to the
environment's python, e.g. a conda env:

```sh
PY=<path-to-env>/bin/python
```

### 3. Run everything

```sh
cd <repo root>
$PY code/run_batch.py                  # matched set (default) -> results/
$PY code/run_batch.py --set legacy     # superseded set        -> results/legacy/
```
Each set is 21 peptides × 3 pairs (63 jobs, about 1–2 minutes), then writes the summary
tables, `RESULTS.md`, and one PyMOL script per peptide. It is resumable: each job writes
`metrics.json` last, and jobs that already have one are skipped. To force a full recompute,
delete the `metrics.json` files. **A set only ever writes inside its own output directory** —
including the PyMOL scripts, which `make_pymol` writes to the morphs directory it is handed —
so running one set cannot overwrite another's results. Job directories are named `<lane>` in
the legacy set and `<lane>_to_<bound>` in the matched set, so the two can share a tree safely. For long runs, `code/run_batch.sh` launches it
detached with `caffeinate`.

### 4. Outputs (`results/`)

```
RESULTS.md                one row per peptide: strain per lane + Apo/Bound patch columns
strain_summary.csv        per (peptide, lane): E_apo, E_boltz, dE, heavy-atom RMSDs, note flags
patch_summary.csv         per (peptide, lane): exposed +/−/hydrophobic counts, buried / newly exposed
skipped.csv               peptides with no co-fold found (currently none)
run.log                   progress log
morphs/                   all morph movies <mol>_<peptide>_<lane>_to_boltz.pdb + one .pml per peptide
<mol>/<peptide>/_boltz/   bound structure relaxed once and shared by all three lanes
<mol>/<peptide>/<lane>/   that lane's job directory (contents below)
```

Per-lane job directory (`<mol>/<peptide>/<lane>/`), one peptide–apo-model pair:

| File | What it is | From |
|---|---|---|
| `morph.pdb` | The apo → bound movie: one multi-MODEL PDB, each MODEL a frame (copied into `morphs/` as `<mol>_<peptide>_<lane>_to_boltz.pdb`) | `morph.py` |
| `<peptide>_<lane>.pdb` | The **apo** structure, patch data ready for coloring: occupancy (`spectrum q`) = residue charge (−1/+1), B-factor (`spectrum b`) = Kyte–Doolittle hydrophobicity (0–1). These are **not colored on load** — a loaded PDB shows as a flat single-color cartoon until a `spectrum` command runs (see below) | `surface_patches.py` |
| `<peptide>_boltz.pdb` | The **bound** co-fold, same data as the apo file, superposed onto the apo frame so the pair overlays in PyMOL | `surface_patches.py` |
| `<peptide>_patches.pml` | The viewer script that does the coloring: load it from this folder (`pymol <peptide>_patches.pml`) and it loads both PDBs, applies the spectrum, and shows the apo state — click the other name in the object panel to switch. At the prompt, `hydro` switches the palette to hydrophobicity, `charge` back to charge. **Charge palette (blue–white–red): blue = negatively charged residue, white = neutral, red = positively charged. Hydrophobicity palette (white–orange): white = polar, orange = hydrophobic.** | `make_pymol.py` |
| `<lane>.pdb` | The apo fold reduced to the peptide chain, as handed to the strain calculation (before hydrogens/relaxation) | `strain.py` |
| `<lane>_H.pdb` | That structure with PDBFixer hydrogens (pH 7) — cached so repeated strain runs are bit-reproducible | `strain.py` |

Plus, in `<mol>/<peptide>/_boltz/`: `boltz.pdb` (the co-fold reduced to the peptide chain)
and `boltz_H.pdb` (its cached protonation) — used by all three lanes so the bound energy is
identical across them. Alongside the PDBs, the lane directory has
`<peptide>_residues.csv` (per-residue charge, KD, rSASA in each state),
`<peptide>_summary.txt` (human-readable patches + changes) and `metrics.json` (all numbers;
written last = job's done-sentinel).

Text outputs (`RESULTS.md`, CSVs, per-job `*_summary.txt` / `*_residues.csv`) are tracked in
git. All PDBs, `.pml` scripts and `metrics.json` are regenerated by `run_batch.py`.

If you load one of the two patch PDBs by hand instead of via `*_patches.pml`, color it with
(substituting the object name PyMOL gave it):
```pymol
spectrum q, blue_white_red, <object>, minimum=-1, maximum=1   # charge (occupancy)
spectrum b, white_orange, <object>, minimum=0, maximum=1      # hydrophobicity (B-factor)
```
Each spectrum overwrites the other's atom colors — the last one run wins, so to see the
other view just run its line again. To turn the patch coloring off entirely (back to plain
gray), `color grey, <object>`.

### 5. View the movies in PyMOL

```sh
cd results/morphs
pymol oxybenzone_ox1_orig_f12.pml
```
Each script loads the peptide's three apo → bound movies as three groups in the object panel
(`esmfold`, `of3apo`, `rf3apo`), aligned on the shared bound structure. All three start
**off**; click a name to show that lane. It opens paused on frame 1 (the apo prediction); the
last frame is the bound structure, and the first and last 5 frames are held still. Play with
`mplay`.

Secondary structure is computed **per frame**: each frame is its own object (inside the
group) with its own `dss`, so helices form or melt as the peptide changes shape. The two
endpoints are the real predicted structures; frames in between are linear interpolations, so
their secondary structure is approximate.

Known issues: the GUI pause button doesn't stop playback, so use `mstop` at the PyMOL prompt.
The scripts turn on `suspend_undo` where the PyMOL build supports it; in builds without it,
long looping sessions slowly accumulate undo history from the per-frame show/hide.

To look at the patch maps instead, go to a job's own folder and load the viewer script:
```sh
cd results/oxybenzone/ox2_esm1_f8/of3apo
pymol ox2_esm1_f8_patches.pml
```
It shows the apo state colored by charge — **blue = negative residue, white = neutral,
red = positive** (the bound co-fold is in the object panel, off, already superposed — click
its name to compare). `hydro` / `charge` at the PyMOL prompt switch the palette to
hydrophobicity (**white = polar, orange = hydrophobic**) and back.

### 6. Reading `RESULTS.md`

- **ΔE esm / of3 / rf3** — strain from each apo model, `E(bound) − E(apo)` in kcal/mol.
  Negative: the peptide pays no strain to take the bound shape. Positive: it has to strain.
  Read the **sign when all three lanes agree**, not the magnitude. Large values are mostly
  differences in prediction geometry: the bound structure is the same for all three lanes,
  yet ΔE can swing by over 100 kcal/mol between apo models. ⚠ marks a `note` flag (see
  strain below); don't use those values.
- **Apo (patches)** — per lane, the number of +, −, and hydrophobic patches on the apo
  surface. A patch is 2+ exposed residues of one type within 8 Å (Cβ–Cβ) of each other.
- **Bound (change on folding)** — per lane, what each apo patch became in the bound
  structure: *buried*, *still exposed*, *broken up* (residues still exposed but no longer
  clustered), plus *new* patches that exist only in the bound state.

---

## Repository layout

```
code/morph.py            morph movie between structures (also the shared PDB/CIF reader)
code/strain.py           apo → bound strain energy
code/surface_patches.py  charged / hydrophobic surface patches, apo vs bound
code/make_pymol.py       per-peptide PyMOL loader scripts
code/run_batch.py        runs all of the above for every peptide and lane; writes RESULTS.md
code/run_batch.sh        detached launcher for run_batch.py
data/                    (unused; inputs are read in place from peptidebuilder / boltzgen_local)
results/                 outputs, see above
```

Each script can also be run on its own; the sections below describe them.

---

## code/morph.py — morph movie

Builds a multi-MODEL PDB that linearly interpolates between an ordered list of structures;
each MODEL is one movie frame.

- Reads `.pdb` and Boltz / BoltzGen / OpenFold `.cif`. Peptide only: for `.cif` it keeps the
  polymer `ATOM` records and drops the ligand `HETATM` block.
- Matches atoms across structures by (residue number, atom name), so atom-order differences
  between programs don't matter; only atoms common to every input are morphed.
- Kabsch-superposes each structure onto the first (on Cα by default) so you see the
  conformational change, not whole-body drift.

```sh
R=$MORPH_PEP_ROOT/octinoxate
BOLTZ=$R/boltz/boltz_results_s3_esm2_f4/predictions/s3_esm2_f4/s3_esm2_f4_model_0.cif
$PY code/morph.py -o morph.pdb --frames 30 --dwell 5 $R/folds/esmfold/s3_esm2_f4.pdb $BOLTZ
# several waypoints in one trajectory (A -> B -> C ...)
$PY code/morph.py -o chain.pdb --frames 24 $R/folds/esmfold/s3_esm2_f4.pdb \
    $R/folds/of3apo/s3_esm2_f4.pdb $BOLTZ
```
Options: `--chain` (default A), `--frames` (steps per segment, 30), `--fit CA|backbone|all`,
`--no-fit`, `--dwell` (still frames held at each waypoint). A single morph file can be
viewed with `load morph.pdb` and the state slider or `mplay`.

**Caveat — linear morph.** Interpolation is straight-line in Cartesian space. When the end
states differ a lot (many apo folds here are 10–25 Å Cα RMSD from the bound structure),
intermediate frames can stretch bonds or pass atoms through each other. That is fine for
visualizing the change but not physical. For geometry-preserving intermediates, PyMOL's own
`morph` command (RigiMOL) is slower but better.

## code/strain.py — apo → bound strain

MM single-point energy difference between the apo and bound peptide conformations, ff14SB +
GBn2 implicit solvent (the force field the `openmm` repo uses for proteins).
`dE = E(bound) − E(apo)`; positive means the bound conformation costs the peptide energy
relative to its apo fold. The ligand is not included.

**Protocol (identical for every structure):**
- Add the missing hydrogens (PDBFixer, pH 7).
- Apply an **extremely limited declash relaxation**: high force target (500 kJ/mol/nm) and
  only 30 minimizer iterations. Its only job is to pull steric overlaps in the predicted
  structures apart. Heavy atoms move about 0.1 Å on clean structures and about 0.3 Å on
  clashy ones, far less than the several-Å change that would turn one state into the other.
- Record each state's heavy-atom RMSD from the predicted coordinates. A move past
  `MAX_HEAVY_RMSD` (0.6 Å) is flagged `moved:<state>`; a state whose energy is still positive
  after the relaxation (a clash it couldn't clear) is flagged `clash:<state>`. Flags go in the
  `note` column of `strain_summary.csv` and show as ⚠ in `RESULTS.md`.
- In the batch, the bound structure is relaxed once per peptide and reused for all three
  lanes, so differences between lanes come only from the apo models.
- Reproducible: PDBFixer places hydrogens slightly differently on each call, so the
  protonated structure is cached (`<label>_H.pdb`) and OpenMM runs single-threaded. Same
  inputs give the same dE.

```sh
$PY code/strain.py --labels esmfold,boltz --workdir /tmp/strain_demo \
    $R/folds/esmfold/s3_esm2_f4.pdb $BOLTZ
```
Options: `--iters` / `--tol` (declash budget), `--ph`, `--workdir` (cache peptide +
protonation there).

## code/surface_patches.py — surface patches

Charged and hydrophobic patches on the apo and bound surfaces. Uses static per-residue charge
(Arg/Lys +1, Asp/Glu −1; His neutral) and Kyte–Doolittle hydrophobicity (≥ 1.5 counts as
hydrophobic), gated by relative solvent exposure (rSASA > 0.25, via mdtraj). Exposed residues
of one type whose Cβ atoms are within 8 Å are grouped into patches, and the two states are
compared.

```sh
$PY code/surface_patches.py --tag s3_esm2_f4 --labels esmfold,boltz --outdir patches \
    $R/folds/esmfold/s3_esm2_f4.pdb $BOLTZ
```
Outputs in `--outdir`:
- `<tag>_<state>.pdb` — peptide with occupancy = residue charge (−1/+1), B-factor =
  hydrophobicity (`spectrum q, blue_white_red, minimum=-1, maximum=1` for charge — blue is
  negative, white neutral, red positive; `spectrum b, white_orange, minimum=0, maximum=1`
  for hydrophobicity — white is polar, orange hydrophobic). On its own this PDB
  loads flat; the spectrum command (or the generated `<tag>_patches.pml`) is what colors it.
- `<tag>_patches.pml` — PyMOL loader written next to the PDBs: applies the spectra and
  shows the pair, with `hydro` / `charge` to switch palettes
- `<tag>_residues.csv` — per residue: charge, KD, rSASA in each state, exposure change
- `<tag>_summary.txt` — patches in each state and what changed

Why static charges rather than per-atom AM1-BCC (OpenFF): AM1-BCC on a whole 34-residue
peptide (~626 atoms) fails (`sqm` errors out, and OpenFF warns it isn't meant for molecules
over ~150 atoms). Static residue charges are what a charged/hydrophobic patch map needs anyway.

## code/make_pymol.py — PyMOL scripts

Writes `results/morphs/<mol>_<peptide>.pml` for every peptide with morphs in
`results/morphs/` (behavior described in Workflow step 5). `run_batch.py` calls it at the end;
to regenerate the scripts alone:
```sh
$PY code/make_pymol.py
```

## code/run_batch.py — everything, every peptide

Finds peptides as the apo folds that have a bound structure (Boltz, else BoltzGen), and for
each peptide × lane runs the morph, strain and patch analysis, then rebuilds the CSVs,
`RESULTS.md` and the PyMOL scripts from all `metrics.json` files present. See Workflow above
for running it and for the output layout.
