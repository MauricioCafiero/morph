# Results — apo → Boltz co-fold strain and surface-patch changes

One row per peptide. **Strain ΔE** = E(co-fold) − E(apo) in kcal/mol, ff14SB + GBn2,
each predicted structure scored as-is with only a limited declash relaxation (heavy
atoms move <0.6 Å; see README). Negative ⇒ the bound conformation is *lower* energy
than the apo fold; positive ⇒ the apo fold must strain to reach the bound state.
ΔE is given for each independent apo model (ESMFold / OF3-apo / RF3-apo).

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
structure is the same for all three lanes; the changes differ because the apo models do.

⚠ = a `note` flag in `strain_summary.csv` (clash or heavy-atom move past guard).

| Mol | Peptide | Kind | ΔE esm | ΔE of3 | ΔE rf3 | Apo (patches) | Bound (change on folding) |
|-----|---------|------|-------:|-------:|-------:|---------------|---------------------------|
| octinoxate | bg33_1 | boltzgen | -128 | -242⚠ | +21 | **esm** 0 +, 0 −, 1 hydrophobic<br>**of3** 0 +, 0 −, 1 hydrophobic<br>**rf3** 0 +, 0 −, 1 hydrophobic | **esm** hydrophobic: 1 still exposed<br>**of3** hydrophobic: 1 still exposed<br>**rf3** hydrophobic: 1 still exposed |
| octinoxate | bg33_2 | boltzgen | -272⚠ | -23⚠ | +4 | **esm** 0 +, 0 −, 1 hydrophobic<br>**of3** 0 +, 0 −, 1 hydrophobic<br>**rf3** 0 +, 0 −, 1 hydrophobic | **esm** hydrophobic: 1 still exposed<br>**of3** hydrophobic: 1 still exposed<br>**rf3** hydrophobic: 1 still exposed |
| octinoxate | bg33_3 | boltzgen | +45 | +59 | +24 | **esm** 1 +, 1 −, 1 hydrophobic<br>**of3** 1 +, 1 −, 1 hydrophobic<br>**rf3** 1 +, 1 −, 1 hydrophobic | **esm** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed |
| octinoxate | bg33_4 | boltzgen | -32 | -48 | -70 | **esm** 0 +, 1 −, 1 hydrophobic<br>**of3** 0 +, 1 −, 1 hydrophobic<br>**rf3** 0 +, 1 −, 1 hydrophobic | **esm** −: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** −: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** −: 1 still exposed; hydrophobic: 1 still exposed |
| octinoxate | orig_f12 | designed | +18 | -66 | -139 | **esm** 0 +, 0 −, 0 hydrophobic<br>**of3** 0 +, 0 −, 0 hydrophobic<br>**rf3** 0 +, 1 −, 0 hydrophobic | **esm** −: 1 new<br>**of3** −: 1 new<br>**rf3** −: 1 broken up, 1 new |
| octinoxate | s2_esm2_control | esm2-variant | +95 | +57 | 0 | **esm** 1 +, 1 −, 1 hydrophobic<br>**of3** 1 +, 1 −, 1 hydrophobic<br>**rf3** 1 +, 1 −, 1 hydrophobic | **esm** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed |
| octinoxate | s3_esm2_f4 | esm2-variant | +1 | -24 | +45 | **esm** 1 +, 0 −, 2 hydrophobic<br>**of3** 1 +, 0 −, 2 hydrophobic<br>**rf3** 1 +, 0 −, 2 hydrophobic | **esm** +: 1 still exposed; hydrophobic: 2 still exposed<br>**of3** +: 1 still exposed; hydrophobic: 2 still exposed<br>**rf3** +: 1 still exposed; hydrophobic: 1 still exposed, 1 broken up, 1 new |
| octinoxate | s3_orig_f12 | designed | -21 | -55 | +18 | **esm** 0 +, 0 −, 0 hydrophobic<br>**of3** 1 +, 0 −, 1 hydrophobic<br>**rf3** 1 +, 0 −, 0 hydrophobic | **esm** no patches<br>**of3** +: 1 broken up; hydrophobic: 1 broken up<br>**rf3** +: 1 broken up |
| octinoxate | shuffle_control | null | +133 | +96 | +69 | **esm** 0 +, 0 −, 0 hydrophobic<br>**of3** 0 +, 1 −, 0 hydrophobic<br>**rf3** 0 +, 0 −, 1 hydrophobic | **esm** −: 1 new; hydrophobic: 1 new<br>**of3** −: 1 still exposed; hydrophobic: 1 new<br>**rf3** −: 1 new; hydrophobic: 1 still exposed |
| octinoxate | shuffle_control_esm0 | null-esm | -103 | -8 | -49 | **esm** 2 +, 3 −, 3 hydrophobic<br>**of3** 3 +, 2 −, 3 hydrophobic<br>**rf3** 3 +, 2 −, 3 hydrophobic | **esm** +: 2 still exposed, 1 new; −: 3 still exposed; hydrophobic: 3 still exposed<br>**of3** +: 3 still exposed; −: 2 still exposed; hydrophobic: 3 still exposed<br>**rf3** +: 3 still exposed; −: 2 still exposed; hydrophobic: 3 still exposed |
| oxybenzone | bgox31_2 | boltzgen | +27 | +55 | +74 | **esm** 1 +, 0 −, 1 hydrophobic<br>**of3** 1 +, 0 −, 1 hydrophobic<br>**rf3** 1 +, 0 −, 1 hydrophobic | **esm** +: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; hydrophobic: 1 still exposed<br>**rf3** +: 1 still exposed; hydrophobic: 1 still exposed |
| oxybenzone | bgox31_3 | boltzgen | -7 | +48 | +39 | **esm** 2 +, 3 −, 2 hydrophobic<br>**of3** 1 +, 3 −, 2 hydrophobic<br>**rf3** 1 +, 3 −, 1 hydrophobic | **esm** +: 2 still exposed; −: 3 still exposed; hydrophobic: 2 still exposed<br>**of3** +: 1 still exposed, 1 new; −: 3 still exposed; hydrophobic: 2 still exposed, 1 new<br>**rf3** +: 1 still exposed, 1 new; −: 3 still exposed; hydrophobic: 1 still exposed |
| oxybenzone | bgox31_4 | boltzgen | +1 | +82 | -82 | **esm** 1 +, 3 −, 2 hydrophobic<br>**of3** 1 +, 3 −, 0 hydrophobic<br>**rf3** 1 +, 3 −, 0 hydrophobic | **esm** +: 1 still exposed; −: 3 still exposed; hydrophobic: 2 buried<br>**of3** +: 1 still exposed; −: 3 still exposed; hydrophobic: 1 new<br>**rf3** +: 1 still exposed; −: 3 still exposed; hydrophobic: 1 new |
| oxybenzone | bgox31_5 | boltzgen | +27 | +90 | +27 | **esm** 1 +, 1 −, 1 hydrophobic<br>**of3** 1 +, 1 −, 3 hydrophobic<br>**rf3** 1 +, 1 −, 1 hydrophobic | **esm** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 3 still exposed, 1 new<br>**rf3** +: 1 still exposed; −: 1 still exposed; hydrophobic: 1 still exposed |
| oxybenzone | ox1_esm1_f8 | esm2-variant | -214 | -74 | -114 | **esm** 1 +, 1 −, 1 hydrophobic<br>**of3** 1 +, 1 −, 1 hydrophobic<br>**rf3** 1 +, 0 −, 2 hydrophobic | **esm** +: 1 still exposed; −: 1 buried; hydrophobic: 1 still exposed<br>**of3** +: 1 still exposed; −: 1 buried; hydrophobic: 1 still exposed<br>**rf3** +: 1 still exposed; hydrophobic: 2 still exposed, 1 new |
| oxybenzone | ox1_orig_f12 | designed | -12 | -41 | -24 | **esm** 1 +, 0 −, 0 hydrophobic<br>**of3** 0 +, 0 −, 0 hydrophobic<br>**rf3** 1 +, 1 −, 1 hydrophobic | **esm** +: 1 broken up; −: 1 new; hydrophobic: 1 new<br>**of3** −: 1 new; hydrophobic: 1 new<br>**rf3** +: 1 broken up; −: 1 still exposed; hydrophobic: 1 still exposed |
| oxybenzone | ox2_esm1_f8 | esm2-variant | -66 | -21 | -1 | **esm** 0 +, 2 −, 2 hydrophobic<br>**of3** 0 +, 2 −, 3 hydrophobic<br>**rf3** 0 +, 2 −, 3 hydrophobic | **esm** −: 2 still exposed; hydrophobic: 2 still exposed<br>**of3** −: 2 still exposed; hydrophobic: 1 buried, 2 still exposed<br>**rf3** −: 2 still exposed; hydrophobic: 1 buried, 2 still exposed |
| oxybenzone | ox2_orig_f8 | designed | +51 | -19 | -41 | **esm** 0 +, 0 −, 0 hydrophobic<br>**of3** 0 +, 0 −, 1 hydrophobic<br>**rf3** 0 +, 0 −, 0 hydrophobic | **esm** −: 1 new; hydrophobic: 1 new<br>**of3** −: 1 new; hydrophobic: 1 still exposed<br>**rf3** −: 1 new; hydrophobic: 1 new |
| oxybenzone | ox2_shuffle | null | +134 | +51 | +39 | **esm** 0 +, 0 −, 0 hydrophobic<br>**of3** 0 +, 1 −, 0 hydrophobic<br>**rf3** 0 +, 1 −, 0 hydrophobic | **esm** no patches<br>**of3** −: 1 broken up<br>**rf3** −: 1 broken up |
| oxybenzone | ox2_shuffle_esm0 | null-esm | -85 | -54 | -158 | **esm** 1 +, 2 −, 2 hydrophobic<br>**of3** 1 +, 2 −, 3 hydrophobic<br>**rf3** 1 +, 2 −, 2 hydrophobic | **esm** +: 1 still exposed; −: 2 still exposed; hydrophobic: 2 still exposed<br>**of3** +: 1 still exposed; −: 2 still exposed; hydrophobic: 3 still exposed<br>**rf3** +: 1 still exposed; −: 2 still exposed; hydrophobic: 2 still exposed |
| oxybenzone | ox3_orig_f8 | designed | -163 | -154 | +7 | **esm** 0 +, 0 −, 0 hydrophobic<br>**of3** 0 +, 1 −, 1 hydrophobic<br>**rf3** 0 +, 0 −, 1 hydrophobic | **esm** hydrophobic: 1 new<br>**of3** −: 1 broken up; hydrophobic: 1 still exposed<br>**rf3** hydrophobic: 1 still exposed |

---

Per-peptide detail (patch residue lists, per-lane CSV, colored PDBs, morph movie)
lives in `results/<mol>/<peptide>/<lane>/`; all morph trajectories are collected
in `results/morphs/`. Full tables: `strain_summary.csv`, `patch_summary.csv`.
