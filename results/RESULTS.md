# Results — apo → Boltz co-fold strain and surface-patch changes

One row per peptide. **Strain ΔE** = E(co-fold) − E(apo) in kcal/mol, ff14SB + GBn2,
each predicted structure scored as-is with only a limited declash relaxation (heavy
atoms move <0.6 Å; see README). Negative ⇒ the bound conformation is *lower* energy
than the apo fold; positive ⇒ the apo fold must strain to reach the bound state.
ΔE is given for each independent apo model (ESMFold / OF3-apo / RF3-apo).

**Surface patches** are from static residue charge (Arg/Lys +, Asp/Glu −) and
Kyte–Doolittle hydrophobicity gated by solvent exposure. 'Charged (co-fold)' lists the
exposed charged residues in the bound state; 'Hydrophobic Δ on fold' is how many
hydrophobic residues bury (−) or newly expose (+) going apo→co-fold, per lane.

⚠ = a `note` flag in `strain_summary.csv` (clash or heavy-atom move past guard).

| Mol | Peptide | Kind | ΔE esm | ΔE of3 | ΔE rf3 | Charged (co-fold) | Hydrophobic Δ on fold (esm/of3/rf3) |
|-----|---------|------|-------:|-------:|-------:|-------------------|-------------------------------------|
| octinoxate | orig_f12 | designed | +18 | -66 | -138 | 3+ (ARG1,LYS7,LYS17); 4− (ASP4,GLU21,ASP26,GLU31) | −5/+0 / −3/+4 / −2/+7 |
| octinoxate | s2_esm2_control | esm2-variant | +95 | +57 | 0 | 8+ (ARG1,ARG4,LYS5,LYS7,LYS9,LYS17,ARG26,ARG32); 3− (GLU3,GLU6,ASP15) | −1/+0 / −1/+0 / −1/+0 |
| octinoxate | s3_esm2_f4 | esm2-variant | +1 | -24 | +45 | 6+ (LYS8,LYS15,LYS16,LYS17,LYS19,LYS22); 2− (ASP4,GLU27) | −6/+0 / −6/+0 / −5/+6 |
| octinoxate | s3_orig_f12 | designed | -21 | -55 | +18 | 2+ (LYS8,LYS16); 2− (ASP4,GLU27) | −9/+0 / −5/+5 / −8/+8 |
| octinoxate | shuffle_control | null | +133 | +96 | +69 | 2+ (LYS6,LYS27); 2− (GLU21,ASP32) | −9/+0 / −5/+10 / −7/+4 |
| octinoxate | shuffle_control_esm0 | null-esm | -103 | -8 | -49 | 8+ (LYS3,LYS6,ARG9,LYS10,LYS15,LYS19,LYS27,LYS30); 8− (ASP11,ASP14,ASP18,GLU21,GLU23,GLU26,GLU31,ASP32) | −1/+1 / −1/+0 / −0/+0 |
| oxybenzone | ox1_esm1_f8 | esm2-variant | -214 | -74 | -114 | 4+ (ARG3,LYS5,LYS7,ARG16); 3− (GLU11,ASP18,ASP25) | −2/+0 / −2/+0 / −1/+4 |
| oxybenzone | ox1_orig_f12 | designed | -12 | -41 | -24 | 2+ (ARG3,LYS5); 3− (ASP9,GLU11,ASP25) | −3/+0 / −2/+4 / −2/+7 |
| oxybenzone | ox2_esm1_f8 | esm2-variant | -66 | -21 | -1 | 2+ (ARG6,LYS11); 6− (ASP3,ASP15,GLU16,ASP27,ASP28,ASP30) | −4/+1 / −3/+0 / −2/+2 |
| oxybenzone | ox2_orig_f8 | designed | +51 | -19 | -41 | 2+ (ARG6,LYS11); 3− (ASP3,GLU16,ASP27) | −6/+0 / −5/+8 / −5/+7 |
| oxybenzone | ox2_shuffle | null | +134 | +51 | +38 | 2+ (ARG1,LYS27); 3− (GLU3,ASP8,ASP25) | −11/+0 / −7/+6 / −5/+4 |
| oxybenzone | ox2_shuffle_esm0 | null-esm | -85 | -54 | -158 | 6+ (ARG1,LYS12,LYS20,ARG21,LYS23,LYS27); 7− (GLU3,ASP7,ASP8,GLU19,ASP25,GLU28,ASP30) | −1/+1 / −1/+0 / −1/+1 |
| oxybenzone | ox3_orig_f8 | designed | -163 | -154 | +7 | 0+; 3− (ASP1,GLU4,ASP16) | −6/+0 / −4/+4 / −3/+4 |

## Skipped (no Boltz co-fold — BoltzGen arms)

octinoxate/bg33_1, octinoxate/bg33_2, octinoxate/bg33_3, octinoxate/bg33_4, oxybenzone/bgox31_2, oxybenzone/bgox31_3, oxybenzone/bgox31_4, oxybenzone/bgox31_5

---

Per-peptide detail (patch residue lists, per-lane CSV, colored PDBs, morph movie)
lives in `results/<mol>/<peptide>/<lane>/`; all morph trajectories are collected
in `results/morphs/`. Full tables: `strain_summary.csv`, `patch_summary.csv`.
