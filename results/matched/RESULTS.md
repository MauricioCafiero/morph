# Results — matched apo → co-fold strain and surface-patch changes

Each column pairs an apo fold with **its own model's** co-fold: Boltz apo → Boltz
co-fold, OF3 apo → OF3 co-fold, RF3 apo → RF3 co-fold. The apo→bound difference is
therefore the ligand's effect within one model, not a change of model. The legacy
set in `results/` morphed all three apo lanes to the single Boltz co-fold, which
mixed the two; it is kept for comparison. ESMFold has no co-fold, so it is absent
here and present only in the legacy set.

For the Boltz pair, BoltzGen peptides (kind = boltzgen: bg33_*, bgox31_*) keep
BoltzGen's own design co-fold as the bound reference — no Boltz co-fold exists for
those sequences — so that one column is matched for 13 of 21 peptides and
BoltzGen-bound for the other 8.

One row per peptide. **Strain ΔE** = E(co-fold) − E(apo) in kcal/mol, ff14SB + GBn2,
each predicted structure scored as-is with only a limited declash relaxation (heavy
atoms move <0.6 Å; see README). Negative ⇒ the bound conformation is *lower* energy
than the apo fold; positive ⇒ the apo fold must strain to reach the bound state.
ΔE is given for each (apo model → bound model) pair in this set.

**Bound reference.** The Boltz co-fold for designed / variant / null peptides; for
BoltzGen peptides (kind = boltzgen: bg33_*, bgox31_*) it is BoltzGen's own design
co-fold from `boltzgen_local/results/` (matched to the apo folds by sequence). Columns
and files still say "boltz" for both.

**Surface patches.** A patch is 2+ solvent-exposed residues of one type whose Cβ atoms
lie within 8 Å of each other: **+** (Arg/Lys), **−** (Asp/Glu), or **hydrophobic**
(Kyte–Doolittle ≥ 1.5). Isolated exposed residues are not counted. **Apo** gives the
patch counts on each apo model's surface. **Bound** says what happened to them in the
Boltz co-fold: *buried* (most of the patch's residues are no longer exposed), *still
exposed* (the patch persists), *broken up* (residues still exposed but no longer
clustered), and *new* (a patch that only exists in the bound state). The bound
bound structure differs per pair, so both ends of each comparison come from one model.

⚠ = a `note` flag in `strain_summary.csv` (clash or heavy-atom move past guard).

| Mol | Peptide | Kind | ΔE boltz | ΔE of3 | ΔE rf3 | Apo (patches) | Bound (change on folding) |
|-----|---------|------|-------:|-------:|-------:|---------------|---------------------------|
| octinoxate | bg33_1 | boltzgen | +19 | +2 | -17 | **boltz** 0 +, 0 −, 1 hydrophobic<br>**of3** 0 +, 0 −, 1 hydrophobic<br>**rf3** 0 +, 0 −, 1 hydrophobic | **boltz** hydrophobic: 1 still exposed<br>**of3** hydrophobic: 1 still exposed<br>**rf3** hydrophobic: 1 still exposed |
| octinoxate | bg33_2 | boltzgen | +33⚠ | -79⚠ | -16⚠ | **boltz** 0 +, 0 −, 1 hydrophobic<br>**of3** 0 +, 0 −, 1 hydrophobic<br>**rf3** 0 +, 0 −, 1 hydrophobic | **boltz** hydrophobic: 1 still exposed<br>**of3** hydrophobic: 1 still exposed<br>**rf3** hydrophobic: 1 still exposed |
| octinoxate | bg33_3 | boltzgen | +15 | +22 | -48 | **boltz** 1 +, 1 −, 1 hydrophobic<br>**of3** 1 +, 1 −, 1 hydrophobic<br>**rf3** 1 +, 1 −, 1 hydrophobic | **boltz** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed |
| octinoxate | bg33_4 | boltzgen | +30 | +350 | -25 | **boltz** 0 +, 1 −, 1 hydrophobic<br>**of3** 0 +, 1 −, 1 hydrophobic<br>**rf3** 0 +, 1 −, 1 hydrophobic | **boltz** −: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** −: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** −: 1 still exposed; hydrophobic: 1 still exposed |
| octinoxate | orig_f12 | designed | +111 | -22 | -57 | **boltz** 0 +, 0 −, 0 hydrophobic<br>**of3** 0 +, 0 −, 0 hydrophobic<br>**rf3** 0 +, 1 −, 0 hydrophobic | **boltz** −: 1 new<br>**of3** +: 1 new; −: 1 new<br>**rf3** −: 1 still exposed |
| octinoxate | s2_esm2_control | esm2-variant | -58 | +5 | +60 | **boltz** 1 +, 1 −, 1 hydrophobic<br>**of3** 1 +, 1 −, 1 hydrophobic<br>**rf3** 1 +, 1 −, 1 hydrophobic | **boltz** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** +: 1 still exposed; −: 1 broken up; hydrophobic: 1 broken up |
| octinoxate | s3_esm2_f4 | esm2-variant | +18 | +32 | -36 | **boltz** 1 +, 0 −, 2 hydrophobic<br>**of3** 1 +, 0 −, 2 hydrophobic<br>**rf3** 1 +, 0 −, 2 hydrophobic | **boltz** +: 1 still exposed; hydrophobic: 2 still exposed, 1 new<br>**of3** +: 1 still exposed; hydrophobic: 2 still exposed<br>**rf3** +: 1 still exposed; hydrophobic: 2 still exposed |
| octinoxate | s3_orig_f12 | designed | +149 | -37 | -90 | **boltz** 0 +, 1 −, 1 hydrophobic<br>**of3** 1 +, 0 −, 1 hydrophobic<br>**rf3** 1 +, 0 −, 0 hydrophobic | **boltz** −: 1 broken up; hydrophobic: 1 broken up<br>**of3** +: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** +: 1 still exposed; hydrophobic: 1 new |
| octinoxate | shuffle_control | null | +174 | -55 | -26 | **boltz** 0 +, 0 −, 1 hydrophobic<br>**of3** 0 +, 1 −, 0 hydrophobic<br>**rf3** 0 +, 0 −, 1 hydrophobic | **boltz** −: 1 new; hydrophobic: 1 still exposed<br>**of3** −: 1 buried; hydrophobic: 1 new<br>**rf3** hydrophobic: 1 still exposed |
| octinoxate | shuffle_control_esm0 | null-esm | -1783⚠ | +108 | +1892⚠ | **boltz** 3 +, 2 −, 3 hydrophobic<br>**of3** 3 +, 2 −, 3 hydrophobic<br>**rf3** 3 +, 2 −, 3 hydrophobic | **boltz** +: 3 still exposed; −: 2 still exposed; hydrophobic: 3 still exposed<br>**of3** +: 3 still exposed; −: 2 still exposed; hydrophobic: 3 still exposed<br>**rf3** +: 2 still exposed, 1 broken up; −: 1 buried, 1 still exposed; hydrophobic: 1 buried, 2 still exposed |
| oxybenzone | bgox31_2 | boltzgen | -47 | -3 | -13 | **boltz** 1 +, 0 −, 1 hydrophobic<br>**of3** 1 +, 0 −, 1 hydrophobic<br>**rf3** 1 +, 0 −, 1 hydrophobic | **boltz** +: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** +: 1 still exposed; hydrophobic: 1 still exposed |
| oxybenzone | bgox31_3 | boltzgen | +106 | -61 | +5 | **boltz** 2 +, 3 −, 1 hydrophobic<br>**of3** 1 +, 3 −, 2 hydrophobic<br>**rf3** 1 +, 3 −, 1 hydrophobic | **boltz** +: 2 still exposed; −: 3 still exposed; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; −: 3 still exposed; hydrophobic: 2 still exposed<br>**rf3** +: 1 still exposed, 1 new; −: 3 still exposed; hydrophobic: 1 still exposed |
| oxybenzone | bgox31_4 | boltzgen | +80 | +50 | +78 | **boltz** 1 +, 3 −, 0 hydrophobic<br>**of3** 1 +, 3 −, 0 hydrophobic<br>**rf3** 1 +, 3 −, 0 hydrophobic | **boltz** +: 1 still exposed; −: 3 still exposed; hydrophobic: 1 new<br>**of3** +: 1 still exposed; −: 3 still exposed; hydrophobic: 2 new<br>**rf3** +: 1 still exposed; −: 3 still exposed; hydrophobic: 1 new |
| oxybenzone | bgox31_5 | boltzgen | -62 | -1 | -19 | **boltz** 1 +, 1 −, 1 hydrophobic<br>**of3** 1 +, 1 −, 3 hydrophobic<br>**rf3** 1 +, 1 −, 1 hydrophobic | **boltz** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 3 still exposed<br>**rf3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed |
| oxybenzone | ox1_esm1_f8 | esm2-variant | -10 | +38 | -3 | **boltz** 1 +, 1 −, 1 hydrophobic<br>**of3** 1 +, 1 −, 1 hydrophobic<br>**rf3** 1 +, 0 −, 2 hydrophobic | **boltz** +: 1 still exposed; −: 1 buried; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 broken up<br>**rf3** +: 1 still exposed; −: 1 new; hydrophobic: 2 still exposed, 1 new |
| oxybenzone | ox1_orig_f12 | designed | -25 | +21 | +19 | **boltz** 0 +, 1 −, 0 hydrophobic<br>**of3** 0 +, 0 −, 0 hydrophobic<br>**rf3** 1 +, 1 −, 1 hydrophobic | **boltz** −: 1 still exposed; hydrophobic: 1 new<br>**of3** no patches<br>**rf3** +: 1 still exposed; −: 1 broken up; hydrophobic: 1 broken up |
| oxybenzone | ox2_esm1_f8 | esm2-variant | +26 | +25 | +136 | **boltz** 0 +, 2 −, 2 hydrophobic<br>**of3** 0 +, 2 −, 3 hydrophobic<br>**rf3** 0 +, 2 −, 3 hydrophobic | **boltz** −: 2 still exposed; hydrophobic: 2 still exposed<br>**of3** −: 2 still exposed; hydrophobic: 1 buried, 2 still exposed<br>**rf3** −: 2 still exposed; hydrophobic: 1 buried, 2 still exposed |
| oxybenzone | ox2_orig_f8 | designed | -12 | +11 | -35 | **boltz** 0 +, 1 −, 0 hydrophobic<br>**of3** 0 +, 0 −, 1 hydrophobic<br>**rf3** 0 +, 0 −, 0 hydrophobic | **boltz** −: 1 still exposed; hydrophobic: 1 new<br>**of3** −: 1 new; hydrophobic: 1 still exposed<br>**rf3** hydrophobic: 1 new |
| oxybenzone | ox2_shuffle | null | +87 | -58 | -69 | **boltz** 0 +, 0 −, 0 hydrophobic<br>**of3** 0 +, 1 −, 0 hydrophobic<br>**rf3** 0 +, 1 −, 0 hydrophobic | **boltz** no patches<br>**of3** −: 1 still exposed; hydrophobic: 1 new<br>**rf3** −: 1 broken up |
| oxybenzone | ox2_shuffle_esm0 | null-esm | +25 | -58 | +125 | **boltz** 1 +, 2 −, 2 hydrophobic<br>**of3** 1 +, 2 −, 3 hydrophobic<br>**rf3** 1 +, 2 −, 2 hydrophobic | **boltz** +: 1 still exposed; −: 2 still exposed; hydrophobic: 2 still exposed<br>**of3** +: 1 still exposed; −: 2 still exposed; hydrophobic: 3 still exposed<br>**rf3** +: 1 still exposed; −: 2 still exposed; hydrophobic: 2 broken up |
| oxybenzone | ox3_orig_f8 | designed | +70 | -199 | -5 | **boltz** 0 +, 1 −, 1 hydrophobic<br>**of3** 0 +, 1 −, 1 hydrophobic<br>**rf3** 0 +, 0 −, 1 hydrophobic | **boltz** −: 1 broken up; hydrophobic: 1 still exposed<br>**of3** −: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** −: 1 new; hydrophobic: 1 broken up |

---

Per-peptide detail (patch residue lists, per-lane CSV, colored PDBs, morph movie)
lives in `results/<mol>/<peptide>/<lane>/`; all morph trajectories are collected
in `results/morphs/`. Full tables: `strain_summary.csv`, `patch_summary.csv`.
