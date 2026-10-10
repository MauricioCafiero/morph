# Peptide morph trajectories

`morph.py` builds a PyMOL-viewable **multi-MODEL PDB** that interpolates between an
ordered list of peptide structures. Each `MODEL` is one movie frame.

- Reads `.pdb` (apo folds) and Boltz/OpenFold `.cif` (co-folds).
- **Peptide only**: for `.cif` it keeps polymer `ATOM` records and drops the ligand
  `HETATM` block; for `.pdb` it reads chain A `ATOM` records (also no `HETATM`).
- Matches atoms across structures by `(residue number, atom name)`, so order
  differences between programs don't matter; only atoms common to every input are morphed.
- Kabsch-superposes each structure onto the first (on Cα by default) so you see the
  conformational change, not whole-body drift.
- Linear Cartesian interpolation, `--frames` steps per segment.

## Run

```sh
PY=/Users/cafierom/miniforge3/envs/openmm-md/bin/python   # conda env with numpy/openmm/mdtraj
R=~/python_mac/peptidebuilder/runs/octinoxate
BOLTZ=$R/boltz/boltz_results_s3_esm2_f4/predictions/s3_esm2_f4/s3_esm2_f4_model_0.cif

# one apo -> bound morph
$PY morph.py -o morph.pdb --frames 30 --dwell 5 $R/folds/esmfold/s3_esm2_f4.pdb $BOLTZ

# chain several waypoints in one trajectory (A -> B -> C ...)
$PY morph.py -o chain.pdb --frames 24 $R/folds/esmfold/s3_esm2_f4.pdb \
    $R/folds/of3apo/s3_esm2_f4.pdb $BOLTZ
```

Options: `--chain` (default A), `--frames` (steps/segment, 30), `--fit CA|backbone|all`,
`--no-fit`, `--dwell` (still frames held at each waypoint).

## View in PyMOL

```
load morph.pdb, morph
set all_states, off          # show one state at a time
mset 1 x<N>                   # or just press the play ▶ button
# nicer look:
hide everything;  show cartoon;  show sticks
set movie_fps, 15
```
Drag the state slider or hit play to watch the peptide go from the apo prediction into
the Boltz-bound conformation.

## Caveat — linear morph

This is straight-line Cartesian interpolation. When two end states differ a lot (the
apo folds here are 10+ Å Cα RMSD from the co-fold), intermediate frames can stretch bonds
or pass atoms through each other — fine for *visualizing* the change, not physical.
For geometry-preserving intermediates use PyMOL's native `morph` (RigiMOL):
```
load apo.pdb, a
load bound.pdb, b      # must have matching atoms
morph mov, a, b, refinement=3
```
That is slower and needs identical atom sets already loaded; `morph.py` is the quick,
dependency-light path and handles cif/pdb mixing and ligand stripping for you.

## Generated demos (peptide s3_esm2_f4, octinoxate)

- `morph_s3_esm2_f4_esmfold_to_boltz.pdb`
- `morph_s3_esm2_f4_of3apo_to_boltz.pdb`
- `morph_s3_esm2_f4_rf3apo_to_boltz.pdb`

---

# strain.py — apo -> cofold conformational strain

MM single-point energy difference between the apo and co-fold peptide conformations,
`ff14SB + GBn2` implicit solvent (the openmm repo's protein force field). `dE = E(cofold)
- E(apo)`; positive means the bound conformation costs the peptide internal energy relative
to its apo fold.

**Protocol (fixed, applied identically to every structure):**
- Score the predicted structure essentially *as-is*. We add the missing **hydrogens**
  (PDBFixer, pH 7), then apply an **extremely limited declash relaxation** whose only job is
  to pull steric overlaps apart: **high force target (500 kJ/mol/nm) + just 30 iterations**.
  Clash gradients (r^-12) are enormous, so the first steps move exactly the overlapping
  atoms; heavy atoms shift ~0.1 Å on clean structures, ~0.3 Å on clashy ones — local
  declashing, not the several-Å move that would turn one state into the other.
- Each state's **heavy-atom RMSD from the predicted coordinates** is recorded; anything past
  `MAX_HEAVY_RMSD` (0.6 Å) is flagged `moved:<state>` in `strain_summary.csv`'s `note` column.
- The **co-fold is relaxed once per peptide** and reused across its three apo lanes, so dE
  differences between lanes come only from the apo models.
- **Reproducible**: the only run-to-run noise was PDBFixer placing hydrogens differently each
  call, so the protonation is cached (`<label>_H.pdb`) and the CPU runs single-thread.

```sh
PY=/Users/cafierom/miniforge3/envs/openmm-md/bin/python     # needs openmm + pdbfixer
R=~/python_mac/peptidebuilder/runs/octinoxate
BOLTZ=$R/boltz/boltz_results_s3_esm2_f4/predictions/s3_esm2_f4/s3_esm2_f4_model_0.cif
$PY strain.py --labels esmfold,boltz --workdir /tmp/strain_demo \
  $R/folds/esmfold/s3_esm2_f4.pdb $BOLTZ
```
`--iters` / `--tol` adjust the declash budget; `--workdir` caches the peptide + protonation.

> The limited declash was added because, with heavy atoms fully frozen, a steric clash
> already present in a prediction cannot relax and blows the energy up (r^-12); several raw
> Boltz co-folds and apo predictions had such clashes. The bounded relaxation removes them
> while keeping the conformation (RMSD < 0.6 Å), so dE reflects the fold, not a contact.

# surface_patches.py — charged / hydrophobic surface spots, apo vs cofold

Where are the charged and hydrophobic patches on the apo surface, and do they move, bury,
or newly form in the co-fold? Uses **static** per-residue charge (Arg/Lys +1, Asp/Glu -1)
and Kyte-Doolittle hydrophobicity, gated by solvent exposure (rSASA via mdtraj). No QM.
Spatially adjacent exposed spots of the same type are grouped into patches and the two
states are diffed.

```sh
PY=/Users/cafierom/miniforge3/envs/openmm-md/bin/python     # needs mdtraj
$PY surface_patches.py --tag s3_esm2_f4 --labels apo_esm,boltz --outdir patches \
    $R/folds/esmfold/s3_esm2_f4.pdb $BOLTZ
```
Outputs in `patches/`:
- `<tag>_<state>.pdb` — peptide with **B-factor = residue charge**, **occupancy = hydrophobicity**
- `<tag>_residues.csv` — per-residue charge / KD / rSASA(apo,cof) / exposure change
- `<tag>_summary.txt` — patches per state + the apo->cofold diff

View the surface in PyMOL:
```
load patches/s3_esm2_f4_boltz.pdb
show surface
spectrum b, blue_white_red, minimum=-1, maximum=1      # charge: blue +, red -
# or hydrophobicity:
spectrum q, white_orange, minimum=0, maximum=1
```

> **Note — why not AM1-BCC charges?** Per-atom AM1-BCC (OpenFF) on the whole 34-mer
> (~626 atoms, net +4) fails: `sqm` fatal-errors and OpenFF itself warns AM1-BCC is not
> meant for >150-atom molecules. Conformation-dependent per-atom charge values are not
> obtainable on the intact peptide, so the surface analysis uses static residue charges —
> which is what "charged vs hydrophobic patches" actually needs.

---

# run_batch.py — all peptides, all lanes

Runs morph + strain + surface patches for every peptide that has a Boltz co-fold, across
the three apo lanes (ESMFold, OF3-apo, RF3-apo). Resumable (a `metrics.json` per job is the
done-sentinel) and reap-proof; launch it detached via `code/run_batch.sh` for long runs.

```sh
/Users/cafierom/miniforge3/envs/openmm-md/bin/python code/run_batch.py
```

Layout under `results/`:
```
<mol>/<peptide>/<lane>/   morph.pdb, <pep>_<lane>.pdb + <pep>_boltz.pdb (patch-colored),
                          <pep>_residues.csv, <pep>_summary.txt, <lane>_H.pdb, metrics.json
<mol>/<peptide>/_boltz/   boltz.pdb, boltz_H.pdb  (co-fold relaxed once, shared by all lanes)
morphs/                   every morph trajectory, named <mol>_<peptide>_<lane>_to_boltz.pdb
RESULTS.md                one-row-per-peptide table: per-lane strain dE + patch changes
strain_summary.csv        per (peptide,lane): E_apo, E_boltz, dE, heavy RMSDs, note flags
patch_summary.csv         per (peptide,lane): exposed +/−/hydrophobic counts, buried/newly exposed
skipped.csv               peptides with no Boltz co-fold (BoltzGen bg33_*/bgox31_*)
```

**Where are the PDB files?** The viewable structures live in two places: the ready-to-load
**morph movies are all in `results/morphs/`** (one file per peptide×lane), and the
patch-colored endpoint structures + per-job movie sit in each
`results/<mol>/<peptide>/<lane>/` folder alongside that job's numbers.
