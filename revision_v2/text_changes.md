# Phase D — text changes, v1 (frozen, `manuscript/*.md`) → v2 (`manuscript/revised_v2/*.md`)

Every number in the new text is a key of `manuscript/revised_v2/generated/numbers.tsv` (written by `revision_v2/scripts/results_numbers.py`; each row names its evidence file). Section-level entries below give: what was removed, what was added or changed, the numbers.tsv keys used, and the reason (reviewer comment or correction). The v1 → v2 diff itself is `diff manuscript/<section>.md manuscript/revised_v2/<section>.md`.

Conventions kept from v1: British spelling, "SNP-mediated reference-genome mapping bias", "naive / WASP / FVB arm", `[CITE Author Year]` placeholders (renumbering in Phase F).

## Abstract
- Removed: five datasets / mammary tumour as a comparison tissue; 32 genes (0.14% of 22,548); old SNP-density medians; the old two worked examples; "cardiac" claims.
- Added/changed: four tissues from two laboratories, sequence-verified samples; 36 FVB/NJ + 34 C57BL/6J animals (`n_animals_strain_comparison_FVB/_B6`); WASP attrition 0.2–2.6% (`wasp_removed_min_pct`, `wasp_removed_max_pct`); ≤ 0.014 (`gw_wasp_minus_naive_absmax`); median within 0.007 (`gw_median_bias_absmax`); catalog 138 = 0.50% of 27,732 (`catalog_n`, `catalog_pct_of_tested`, `genes_tested_union`), 46 in all four (`catalog_n_tissues_4`); SNP density 7.19 vs 0.25, p = 5.2 × 10⁻⁵⁶ (`pred_catalog_median_snp_density`, `pred_noncatalog_median_snp_density`, `pred_catalog_vs_noncatalog_snp_density_MWU_p`); Ndufs5 3.3–10.2 (`ex_robust_*_bias`); Rps2 59-fold / 1.5-fold (`ex_extreme_naive_fold`, `ex_extreme_fvb_fold`); negative control and within-strain sentences. F1 numbers unchanged (re-derived: `f1_naive_pooled_mean/_se/_p`).

## Introduction
- Paragraph 1: added one sentence on the scale of FVB/NJ–B6 divergence (R1-1): ~5.0 million homozygous autosomal SNPs, 1.2 million indels, 203,510 exonic SNPs (`vcf_fvbnj_homalt_autosomal_snps`, `vcf_fvbnj_indels`, `vcf_fvbnj_homalt_exonic_snps`; evidence `revision_v2/audit/fvb_vcf_counts.txt`), and the non-uniform distribution (gene-level densities in `gene_meta.tsv`).
- Paragraph 2: added the clause on the published FVB/NJ assembly [CITE Wong 2012] (R1-2/R1-3). Paragraph 4: added the sentence on sequence verification of public samples. Paragraph 5 (summary): four tissues / two laboratories / 36 + 34 animals / negative control; 0.2–2.6%; 138 genes; Ndufs5 and Rps2 replace the old examples; 59-fold (`ex_extreme_naive_fold`). Removed every reference to the old tissue set and to the old examples.

## Results
- §1 Survey: unchanged except "(Fig. 1B)" added (R1-4).
- §2 F1 ASE: dataset description per manifest (`n_f1`; GSE293440 attribute "B6xFVB", cerebral cortex); "(Fig. 2A–C)" cites all panels; sentence that all 12 animals passed the identity check. Numbers re-derived and identical to v1 (`f1_naive_pooled_mean` 0.3329, `f1_naive_t` −43.70, `f1_naive_p` 1.1 × 10⁻¹³, `f1_wasp_pooled_mean` 0.3856, `f1_wasp_t` −12.33, `f1_wasp_p` 8.8 × 10⁻⁸, `f1_naive_sites_min/max/mean`, per-dataset means).
- §3 NEW "Every sample entering the pure-strain comparison was verified by sequence": manifest 186 rows (`n_manifest_rows`); 83 animals matched (`qc_used_rows` − `n_f1`; `qc_verdict_B6_pure` 47, `qc_verdict_FVBNJ_pure` 36; `qc_B6_mean_af_max` 0.018; `qc_FVB_mean_af_min` 0.987; `qc_pure_n_sites_min/max` 112–3,531); F1 (`qc_F1_mean_af_min/max` 0.474–0.481, `qc_F1_frac_mid_min` 0.955); GSE55792 rejection (`qc_failed_mean_af_min/max` 0.181–0.203, `qc_failed_frac_hi_min/max` 14–18%) with Supplementary Fig. S1 and neutral wording ("fixed non-B6 haplotype blocks … consistent with an incompletely backcrossed background"; decision 2026-10-14). Rule 0.4-4.
- §4 Three-arm comparison: tissue list and n per tissue (`n_<tissue>_FVB/_B6`); covariate design; genes tested 27,732 / 18,104–23,047 / 14,316 (`genes_tested_union`, `gw_*_n_genes`, `genes_tested_all4`); medians (`gw_*_median_bias`, `gw_median_bias_absmax`); means (`gw_mean_bias_min/max`); % ≥ 0.5 per tissue (`gw_hypothalamus_pct_ge_0.5` 0.56%, 122; `gw_adipose_pct_ge_0.5` 0.92%, 213); 326 = 1.18% (`genes_ge_0.5_any_tissue`, `_pct`); Fig. 3A representative tissue = liver by rule (`fig3A_tissue`). Fig. 3A now cited before 3B, C (R1-5).
- §5 WASP: 0.23–2.55% (`wasp_removed_min/max_pct`); B6 > FVB in every tissue (`wasp_removed_B6_more_than_FVB_tissues` = 4; `wasp_removed_<tissue>_<strain>_mean_pct`); median SNP density 0.26 (`pred_median_snp_density_all_tested`); Δ per tissue (`gw_*_wasp_minus_naive_mean_logFC`); FVB-arm shift (`gw_fvb_minus_naive_min/max`). Removed the old "0.2–2.3%" and "≤ 0.005" statements.
- §6 Predictors: ρ 0.337 (`pred_snp_density_per_kb_rho`; p underflows to 0 → written "p < 10⁻¹⁵"); ρ 0.506 for the 14,316 all-four-tissue genes (`pred_snp_rho_all4tissue_genes`, `pred_n_all4tissue_genes`; decision 2026-10-20 item 1: one sentence); 117/138, q75 2.56 (`pred_catalog_in_top_quartile`, `pred_q75_snp_density`); top decile ≥ 6.79 (`pred_q90_snp_density`; `pred_top_vs_bottom_decile_bias_MWU_p` underflow); length ρ 0.051, p 2.7 × 10⁻¹⁷; expression ρ −0.063, p 1.6 × 10⁻²⁵; binned mean within ±0.07 (`pred_expression_binned_mean_absmax`); biotype means (`pred_biotype_*_mean`), medians within ±0.003 (`pred_biotype_median_absmax`); pseudogene interpretation reworded (observation (e) of the Phase C summary).
- §7 Catalog: 138 (0.50%); control-only 132 / 124 shared (`catalog_water_n`, `catalog_water_shared`); 46/35/57 (`catalog_n_tissues_*`); per-tissue counts (`catalog_biased_in_*`); biotypes 102/15/18/3 (`catalog_biotype_*`, `catalog_other_biotypes`); SNP density 7.19 vs 0.25, p 5.2 × 10⁻⁵⁶; top decile 72/138, OR 10.0 (7.2–14.1), p 1.1 × 10⁻³⁵ (`pred_top_decile_*`); 16 genes with mean < 0.5 (`catalog_mean_bias_lt_0.5`; observation (a)); Ndufs5 (`ex_robust_*`); named four-tissue genes read from `catalog/catalog.primary.tsv` (n_tissues = 4); Gm25939 (`catalog.primary.tsv` row). Fig. 5A now cited. Removed: all v1 gene names and values, the substrain-support columns (every animal is now FVB/NJ), the PyVT-discordance paragraph.
- §8 NEW negative control: PyVT n = 13 (`n_pyvt`), mean AF 0.004–0.018 (Supplementary Table S1 rows); genome-wide median 0.000 and 0.3–0.5% (`nc_*_gw_median_mean`, `nc_*_gw_pct_ge_0.5_mean`); catalog medians FVB +0.24 to +0.49, B6 −0.39 to −0.69, PyVT −0.08 to −0.16 (`nc_<group>_catalog_median_min/max`), 77 genes (`nc_pyvt_catalog_genes_n`); held-out B6 −0.43 to −0.66, FVB +0.29 to +0.44, PyVT −0.09 to −0.18 (`ho_*_median_min/max`); chr 16 block 37.7–49.7 Mb in 13/13 and 0/16 (`fs_chr16_block_*`); Cd200r4 −0.26, Gm6030 not expressed, median without them −0.12 (`nc_pyvt_median_swap_Cd200r4`, `nc_pyvt_median_swap_Gm6030`, `nc_pyvt_catalog_median_excl_chr16block_mean`). Decision 2026-10-20 item 2 (wording written after the results, as required).
- §9 Multi-mapping: 82,264 pairs, median +0.008, mean +0.021, 1.1% > 0.5, 0.7% > 1 (`mm_*`); 82/56 classes, mean |Δ| 1.01/0.17, 37 ≤ 0.30, 56-class residual 0.36–1.14, 19 ≥ 0.80 (`mech_*`). Panels S2A–D cited.
- §10 Worked examples: Ndufs5 and Rps2 chosen by rule (`ex_robust_*`, `ex_extreme_*`; `worked_examples.tsv`); all per-tissue values from `ex_<role>_<tissue>_*`; Rps2 skin fold 59 (`ex_extreme_naive_fold`), 1.5 (`ex_extreme_fvb_fold`), per-animal CPM means (`ex_extreme_cpm_log2_*`). Fig. 6B now cited before 6D. Decision item 5: one gene per mechanism class.
- §11 NEW within-strain comparisons (R2-1): `wf_*` keys (DE counts, concordance, r, % |Δ| ≥ 0.5, catalog median |Δ|).

## Discussion
- Deleted: the adipose mapping-rate paragraph (based on misidentified samples) and the PyVT CNV / pvt_discordant paragraph (the dataset is a B6 congenic; [CITE Lin 2003] dropped).
- Merged the three WASP paragraphs into one (R2-3; keeps the F1 vs pure-strain contrast, the sparsity argument, and "not a failure of WASP"), with the new numbers and the Fig. 6B observation on how WASP acts at catalog loci.
- Substrain paragraph rewritten: every pure-strain animal is now FVB/NJ (JAX); the F1 datasets' substrain status.
- Mechanism paragraph: Ndufs5 / Rps2 and the 82/56 counts (`mech_*`); practical implication (37 of 82).
- NEW paragraph on the FVB/NJ assembly question (R1-2/R1-3): why a SNP-patched GRCm39 was used, when FVB_NJ_v3 is preferable, and that the direct comparison (optional Phase G) was not performed.
- NEW paragraph on sequence-based verification of public samples (GSE55792 example with neutral wording; PyVT chr 16 segment; recommendation) — Appendix D of the instructions.
- NEW paragraph on the negative control, including the unexplained smaller magnitude (`nc_PyVT_SE_catalog_median_min/max` vs `nc_C57BL6J_*`), stated without interpretation.
- NEW paragraph on within-strain comparisons (R2-1).
- NEW "Conclusions" section with the substrain caveat (R2-2) and the identity-check recommendation.
- Numbers: 138 / 0.50% / 27,732; 0.007; −0.02 to −0.03 (`gw_mean_bias_min/max` rounded); 0.23–2.55%; 0.26; 0.014; 1.19–2.55% vs 0.23–0.37% (`wasp_removed_<tissue>_<strain>_mean_pct`); 7.19 vs 0.25; residual fractions (`ex_extreme_*_residual`, `ex_robust_residual`).

## Methods
- "RNA-seq datasets": replaced entirely by `generated/methods_datasets.md` (rendered from the manifest; marked with BEGIN/END comments). The v1 sentences "confirmed from GEO series metadata" / "recorded in GEO sample attributes" are gone (Appendix C).
- NEW "Sample identity verification" (site set 203,510 / subset 20,351 with seed; ASEReadCounter parameters; verdict rules from `analysis_v2/pipeline/identity_verdict.py`; validation of the subset from `revision_v2/speedup/03_qc_subsample.md`; Supplementary Table S1 / Fig. S1).
- "Personalised reference": gene counts 1,727 of 56,941, all with zero SNPs (`gtf_genes_*`); STAR index overhang stated.
- "Three-arm alignment", Arm 3: B6 GTF inserted at mapping time; effective junction set (`junctions_*`; Supplementary Note 1); scheduler sentence; Supplementary Table S4.
- "Quantification": paired-end datasets listed per manifest; single-end GSE86315 and GSE175625.
- "Differential expression": designs `~ diet + strain` and `~ treatment + sex + strain` (the sample-sheet column is called `diet` for both studies in the code); control-only sensitivity with n per tissue (from `generated/table1.md`).
- "Catalog construction": four tissues; the fixed-criterion statement; matched-reference sensitivity (decision item 4) → Supplementary Table S3c; removed the PyVT exclusion / pvt_discordant / FVB-N support columns and the v1 example gene.
- NEW "Negative control" (reference-swap ratio, held-out comparison, chr 3/10/16 full-site scan with comparators).
- "Predictor analysis": definition = mean over tested tissues, every gene tested in ≥ 1 tissue (decision item 1), plus the all-four-tissue correlation.
- "Multi-mapping": per tissue, same models.
- NEW "Choice of worked examples and within-strain comparison".
- "Statistical analysis": tests listed; PDF + 600-dpi TIFF; data behind every panel.
- "Availability": GSE86315 added, GSE135230 removed; `v2-corrected` branch/tag; manifest, QC and change log released.

## Figure legends
- All legends rewritten for the v2 figures (`analysis_v2/figures/*_data.tsv` are the sources: Fig1_survey, Fig2_f1_ase, Fig3_strain_comparison, Fig4_predictors, Fig5_catalog, Fig6_worked_examples, SFig1_GSE55792_haplotype_blocks, SFig2_multimapping). Supplementary figures renumbered: S1 = identity check of GSE55792 (cited first), S2 = multi-mapping.

## Table 1
- Replaced by `generated/table1.md` (make_table1.py; QC column now filled: 14/14, 16/16, 24/24, 16/16; 13/13; 8/8, 4/4).

## Supplementary tables and notes
- Captions rewritten for S1 (manifest + identity verification), S2 (catalog, with `in_control_only_catalog`, `in_matched_reference_catalog`, `pyvt_median_reference_swap`), S3 (a–f: sensitivity grids incl. matched-reference, negative control per animal / held-out / PyVT segments), S4 (alignment stats per animal × arm + attrition by group), S5 (survey; was S1). Tables built by `revision_v2/scripts/build_supplementary_tables.py` into `manuscript/revised_v2/supplementary_tables/`.
- Supplementary Note 1: v1 validation table (old samples, old catalog) removed; added the effective junction set section (`revision_v2/audit/junction_sets_fvb_index_vs_b6_gtf.txt`) and the B6-GTF-at-both-stages description.
- Supplementary Note 2: all numbers and genes replaced (`mm_*`, `mech_*`, `ex_*`); the Rps2 skin case replaces the v1 example.

## References
- [CITE Lin 2003] removed (CNV paragraph deleted); [CITE Wong 2012] added (DOI from the reviewer's comment: 10.1186/gb-2012-13-8-r72). Numbering to be redone in Phase F.

## Appendix C grep — residual hits for human decision (2026-10-20)
- `+1.14` and `0.005` occur in the new Results as genuinely new values (Rps2 liver bias +1.145 → "+1.14", `ex_extreme_liver_bias`; liver WASP − naive mean log₂FC +0.00495 → "+0.005", `gw_liver_wasp_minus_naive_mean_logFC`). They coincide with forbidden v1 strings by chance; left as they are and reported.
- "heart" is the GEO tissue attribute of the rejected GSE55792 samples: it appears in Supplementary Table S1a/S1b (tissue column) and in the row labels of Supplementary Fig. S1 (`SFig1_GSE55792_haplotype_blocks_data.tsv`); the main text says "myocardial". "cardiac" appears only inside GEO series titles in Supplementary Table S5. These data-table hits are reported as EXEMPT by `appendix_c_grep.sh` (not failures); human decision whether to relabel.
