# Strain-identity QC — calibration record (Phase A-2)

Date: 2026-10-13/14. Sites: `analysis_v2/qc/fvbnj_homalt_exonic.hetGT.vcf.gz` (203,510 FVB/NJ homozygous-alt, PASS, biallelic, autosomal, exonic SNPs from MGP REL-2021; GT rewritten to 0/1 only so that GATK ASEReadCounter accepts the sites). Counting: GATK 4.6.2.0 ASEReadCounter, min depth 20, MAPQ ≥ 20, BQ ≥ 20, on the naive (B6-reference STAR) BAMs. Metrics: `revision_v2/scripts/strain_identity_qc.py`. Verdict thresholds were fixed before any sample was scored (instructions v2, A-2) and were not changed afterwards.

| Verdict | Rule |
|---|---|
| FVB/NJ pure | mean_af > 0.90 and frac_hi > 0.85 |
| B6 pure | mean_af < 0.10 and frac_lo > 0.85 |
| F1 (B6 × FVB) | frac_mid > 0.50 and 0.25 ≤ mean_af ≤ 0.55 |
| segregating mixed background | frac_hi, frac_mid, frac_lo all > 0.15 |
| other inbred strain | frac_hi > 0.25 and frac_lo > 0.25 and frac_mid < 0.15 |
| needs review | none of the above |

## Positive controls (identity certain from GEO sample attributes)

| Control | Samples | n_sites (range) | mean_af | frac_hi / frac_mid / frac_lo | Verdict | Expected |
|---|---|---|---|---|---|---|
| B6 (fructose study, GEO strain C57BL/6J) | 18 runs, 6 animals (CFM9/15, CFHyp10/11, CFL10/11) | 4,258–18,436 | 0.001–0.004 | ≤0.001 / ≤0.001 / ≥0.998 | B6 pure ×18 | B6 pure ✔ |
| F1 (GSE200632, GEO "C57BL/6:FVB/N (F1:50/50)") | 4 runs | 13,337–15,207 | 0.474–0.477 | 0.013–0.014 / 0.955–0.958 / 0.029–0.031 | F1 ×4 | F1 ✔ |
| FVB/NJ (fructose study, GEO strain FVB/NJ) | 12 runs, 6 animals | 5,168–18,159 | 0.983–0.993 | ≥0.975 / ≤0.022 / ≤0.006 | FVB/NJ pure ×12 | FVB/NJ pure ✔ |

Both mandated controls (B6, F1) and the FVB/NJ samples fall where the rules expect them; the three classes are separated by > 0.4 in mean_af with no overlap. The old analysis reported an F1 reference-allele fraction of about 0.35 at a different site set (all FVB_NJ SNPs incl. heterozygous/low-quality calls, WASP-filtered BAM); at PASS homozygous exonic sites on the naive BAM, the alt fraction is 0.47–0.48, i.e. a reference bias of ~0.02–0.03 below 0.5. The F1 verdict does not depend on this difference.

## Behaviour on the problem samples

- DBA/2J (GSE123875 DWM6/DWM8): mean_af 0.55–0.60, frac_hi 0.54–0.57, frac_lo 0.36–0.44, frac_mid 0.03–0.07 → "other inbred strain" ×6. Consistent with an inbred strain sharing roughly half of the FVB/NJ-vs-B6 alleles.
- B6.FVB-Tg(MMTV-PyVT) congenic (GSE175625): mean_af 0.011–0.016, frac_lo 0.972–0.982 → "B6 pure" ×13 (the residual 1–2 % is compatible with a small FVB-derived congenic interval around the transgene plus mapping noise).
- GSE135230 Het-cKO (Bin1flox/wt; Thy1-Cre+): mean_af 0.443–0.457, frac_mid 0.925–0.937 → "F1 (B6 × FVB)" ×16. These animals are the first-generation cross of a B6.129S6 floxed line with the FVB/N-Tg(Thy1-cre) line, so they are genuinely F1-like at FVB-informative sites.
- GSE135230 Homo-cKO (Bin1flox/flox; Thy1-Cre+): mean_af 0.237–0.259, frac_mid 0.50–0.56, frac_lo 0.43–0.49, frac_hi ≤ 0.011 → "F1" ×4 and "needs review" ×2. The pattern (half heterozygous, half B6-homozygous, no FVB-homozygous) is that of a backcross (N2) generation, not an F1. The fixed rule set has no class for an N2 and the "segregating mixed background" rule (which needs frac_hi > 0.15) cannot fire. **Reported as is; the rule set was not modified.** Suggested addition for the human to decide (Phase A STOP): "backcross to B6: frac_hi < 0.15 and frac_mid > 0.35 and frac_lo > 0.35".

## Limitation
The check is informative only for FVB/NJ-vs-reference alleles. It distinguishes FVB/NJ, B6, F1 and "not FVB and not B6", but it cannot identify which other strain a sample is (DBA/2J here is inferred from GEO, confirmed only as "non-FVB, non-B6 inbred"), nor FVB/N from FVB/NJ.
