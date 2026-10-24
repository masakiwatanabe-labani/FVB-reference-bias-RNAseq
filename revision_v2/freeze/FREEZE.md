# Number freeze — Phase M (2026-10-24)

Sealed in the commit that adds this file (parent ?, branch v2-corrected); that commit is the one the authors tag v2-corrected-r1.

## What was verified

- **M-0 inputs frozen**: `inputs.sha256` lists every frozen input (numbers.tsv, supplementary_tables_v2.xlsx, every file of RESULTS.md5 and FIGURES.md5, and the two md5 lists); final_check fails if any changes. `build_submission.sh` does not rewrite the xlsx or the figures once frozen.
- **M-1 numeric-token audit** (`token_audit.py`, `resolve_tokens.py` → `token_audit.tsv`): 1784 numeric tokens in the main text (all sections, Table 1, legends), Additional file 1, both letters and the text of the nine figure PDFs.
- **M-2 independent recomputation** (`recompute_independent.py` → `recompute_report.tsv`): no import of the pipeline or of the text/number generators, no numbers.tsv; primary result files only (edgeR tables, count matrices and featureCounts summaries, STAR junction lists and Log.out, F1 and identity per-site tables, WASP allele-removal table, within-strain DE tables, FVB_NJ VCFs and the Ensembl 110 GTF via bcftools/bedtools, rps2_* files, vW tag counts, survey table, submitted-version catalog for the letters). Bounds pass only if the exact value satisfies the bound as printed.
- **D3 double-rounding scan** (`double_rounding_scan.py`): every key printed in the text compared with the value rounded once from full precision (`generated/numbers_fullprecision.tsv`).
- **Seal**: `deliverables.sha256` (every file in revision_v2/submission/, figures included).

## Token classes (M-1)

| class | tokens |
|---|---|
| KEY | 912 |
| STRUCT | 852 |
| M2 | 17 |
| UNMAPPED | 2 |
| XLSX | 1 |

KEY = equals a numbers.tsv value at the printed precision; XLSX = derivable from the supplementary tables; STRUCT = citation, label, accession, year, version, parameter or fixed threshold; M2 = not in numbers.tsv or the tables but reproduced by a PASS claim of the independent recomputation (claim id in `ref`); UNMAPPED = see below.

UNMAPPED (2): **1,047** (GEO search hits, Methods) — carried from the submitted version; no search log exists in the project; not recomputable. Plausibility check 2026-10-23 through E-utilities (db=gds, `FVB AND "expression profiling by high throughput sequencing"[DataSet Type] AND "Mus musculus"[Organism] AND gse[Entry Type]`): 1,073 series today (≥ 1,047, as expected for a later date); text unchanged. **290** (size of the GSE293440 series, Methods) — verified against the GEO SOFT record (290 `!Series_sample_id` entries, 2026-10-23).

## Claims (M-2)

| group | PASS | FAIL | not recomputable |
|---|---|---|---|
| A variant set and references | 21 | 0 | 0 |
| B genome-wide | 24 | 0 | 0 |
| C catalog and predictors | 56 | 0 | 0 |
| D multi-mapping | 31 | 0 | 0 |
| E decomposition | 24 | 0 | 0 |
| F worked examples | 75 | 0 | 0 |
| G F1 allele-specific expression | 43 | 0 | 0 |
| H negative control | 31 | 0 | 0 |
| I WASP attrition | 9 | 0 | 0 |
| J within-strain | 13 | 0 | 0 |
| K letters and survey | 20 | 0 | 2 |
| L Supplementary Fig. S1 | 7 | 0 | 0 |
| N figure text (Phase N) | 17 | 0 | 0 |
| P Phase P (genome-level differences, clusters) | 13 | 0 | 0 |
| Q Phase Q (decomposition wording, within-strain wording) | 15 | 0 | 0 |
| **total** | **399** | **0** | **2** |

## Authorised exceptions to the freeze (authors' decisions, 2026-10-23)

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| D1 | Fig. 1 regenerated (v2 F1 table) | `analysis_v2/figures/Fig1_survey.pdf` | `2c37ea19711eda5b…` | `f71518378c8c2694…` |
| D1 | Fig. 1 regenerated | `analysis_v2/figures/Fig1_survey.tiff` | `c10dc16e6fbc38b9…` | `b144c239491f8251…` |
| D1 | Fig. 1 (submission copy) | `revision_v2/submission/figures/Fig1_survey.pdf` | `2c37ea19711eda5b…` | `f71518378c8c2694…` |
| D1 | Fig. 1 (submission copy) | `revision_v2/submission/figures/Fig1_survey.tiff` | `c10dc16e6fbc38b9…` | `b144c239491f8251…` |
| D1b | Suppl. Fig. S1 regenerated (v2 reference animals) | `analysis_v2/figures/SFig1_GSE55792_haplotype_blocks.pdf` | `e5d6642cc519cead…` | `3374ea196bf338ce…` |
| D1b | Suppl. Fig. S1 regenerated | `analysis_v2/figures/SFig1_GSE55792_haplotype_blocks.tiff` | `c68a08777c22c2b6…` | `a807006a5db1aa29…` |
| D1b | Suppl. Fig. S1 data | `analysis_v2/figures/SFig1_GSE55792_haplotype_blocks_data.tsv` | `8b46b0a4a84629da…` | `6e21a0919d604a90…` |
| D1b | Suppl. Fig. S1 (submission copy) | `revision_v2/submission/figures/SFig1_GSE55792_haplotype_blocks.pdf` | `e5d6642cc519cead…` | `3374ea196bf338ce…` |
| D1b | Suppl. Fig. S1 (submission copy) | `revision_v2/submission/figures/SFig1_GSE55792_haplotype_blocks.tiff` | `c68a08777c22c2b6…` | `a807006a5db1aa29…` |
| D1b | added: reference-animal candidates (+6 per-site tables) | `analysis_v2/results/identity/full_site_tables/candidates.tsv` | `(absent)…` | `0a164012e848db18…` |
| D2 | added: Rps2 exonic variants, union definition | `analysis_v2/results/worked_examples/rps2_exonic_variants_union.tsv` | `(absent)…` | `fbf5d346b6f3cbc5…` |
| D2 | unchanged (no new SNP inside the copied part) | `analysis_v2/results/worked_examples/rps2_pseudogene_alleles.tsv` | `8cf1d351f6a169a0…` | `8cf1d351f6a169a0…` |
| D2/D3/D4/D5/E1 | numbers.tsv (3 values, 59 CPM precisions, new keys) | `manuscript/revised_v2/generated/numbers.tsv` | `b7b3af78ef8afd3f…` | `75b9898203802955…` |
| D1/D1b | FIGURES.md5 (5 lines) | `analysis_v2/figures/FIGURES.md5` | `a15db3c5e3fb67f1…` | `81b7ed8f13cc3428…` |
| D1b/D2 | RESULTS.md5 (8 lines added) | `analysis_v2/results/RESULTS.md5` | `402e3f9f84bfc3d5…` | `128502d394c1f9f1…` |
| D5b | Table 1 tissue cell | `manuscript/revised_v2/table1.md` | `71ab84cee173b941…` | `dd305f1db254fb99…` |
| M-0 | frozen-input list updated for the files above | `revision_v2/freeze/inputs.sha256` | `351616284102d0a4…` | `6c045f5b4c935228…` |

Full before/after hashes: `before_hashes.tsv` (before) and `inputs.sha256` / `deliverables.sha256` (after). Every other file of FIGURES.md5 and RESULTS.md5 is byte-identical to the freeze of 2335065; `rps2_pseudogene_alleles.tsv` is unchanged (no SNP inside the copied part is new).

- **D1** Fig. 1 regenerated (PDF, TIFF; `_data.tsv` byte-identical): `fig1_survey.py` read the submitted-version F1 table; panel B now prints 0.480 / 0.512 (pdftotext), no 0.333 / 0.386.
- **D1b** Supplementary Fig. S1 regenerated: reference animals by rule (most sites of depth ≥ 20 per strain among the v2 animals evaluated at the full site set; candidates in `analysis_v2/results/identity/full_site_tables/candidates.tsv`): adipose_B6_CWM8 (35,058 sites) and adipose_FVB_FWM9 (30,685), both GSE123875; the six per-site tables (3.1 MB in all) are in that directory, the only input of `sfig_gse55792_blocks.py` besides the reference .fai.
- **D2** Rps2 exonic variants recounted with exonic = union of the exons of all Ensembl 110 Rps2 transcripts (`analysis_v2/results/worked_examples/rps2_exonic_variants_union.tsv`): 16 SNPs, 5 indels; 5 SNPs inside the part copied by the pseudogenes (ps10 5/5, ps13 5/5, ps6 3/5 FVB/NJ), 11 outside; 0 indels inside (case 1). The genomic blastn cross-check differs only at 17:24,940,643 (a 106-bp exon below e < 1e-50 as a standalone HSP; inside the copied part in the mature-sequence alignment).
- **D3** Rps2/pseudogene CPM keys stored at full precision; Note 2 table 209.75 → 209.76. Scan: the only double-rounding defect before the change was this value; after: 0.
- **D4 / E2 / M-3a–c / D5 / E1** text changes (trigger PHASE-M; text_changes.md, Phase M). New numbers.tsv keys: `vcf_fvbnj_homalt_autosomal_indels`, `identity_subset_validation_n_animals` (14) and `_max_absdiff` (0.0135), `nc_catalog_genes_expressed_min/max` (76/86), `decomp_sum_mismatch_n_genes` (1) and `decomp_Gm26740_sum_minus_mean_bias` (−0.0186), `sfig1_shared_block_chromosomes` (2, 3, 8, 10 and 11), `sfig1_shared_block_n_chromosomes`, `sfig1_ref_b6_max_hom_fvb` (0.0326), `sfig1_ref_fvb_min_hom_fvb` (0.691). Changed values: `rps2_snps_exonic` 13 → 16, `rps2_indels_exonic` 4 → 5, `rps2_exonic_snps_outside_pseudogene_sequence` 8 → 11; 59 CPM keys changed in stored precision only.
- **D6** `manuscript/revised_v2/generated/wasp_attrition_by_group.tsv` (input of Sheet S4b): input outside analysis_v2/results/, v2-derived, independently verified (M-2 group I, 9/9 PASS); tracked in git; the xlsx was not rebuilt.

## Other records

- Duplicate keys in numbers.tsv (`ho_adipose/hypothalamus/liver/skin_n_catalog`, three rows each; `pred_biotype_kruskal_p`, two rows): identical values in all rows (E3); left unchanged during the freeze.
- D1.1 grep after the changes: no figure or table script reads outside analysis_v2/results/, revision_v2/manifest/ and ref/ except Sheet S4b (D6).
- `revision_v2/freeze/cache/` (intermediate BEDs of M-2, 431 MB) is not committed.
- `submission_checklist.md` was generated before the seal; its seal line therefore reads "not yet written"; the final_check run after the seal (final_check.md) verifies it.

## Phase N (2026-10-23): fixes from the package check — authorised exceptions

No number changed: every _data.tsv is byte-identical to the frozen version; the M-1 token audit of the figure PDFs differs only by layout tokens (axis ticks, digits of rotated labels split by pdftotext), the Fig. 3A key (117 of 138 catalog genes tested in liver; n = 40/31/46 tested in liver) and Fig. 4 p values printed as "p < 1e-15" (was 0.0e+00, 2.7e-17, 1.6e-25, 1.8e-35); all are claims of the M-2 report (group N).

- N-1a Fig. 6D title: "largest single-tissue bias of the direct-misalignment class (skin +5.36 log₂)"; Supplementary Fig. S2D title: "two largest |Δ| of the catalog" (the only other "largest" in a figure).
- N-1b Fig. 4: printed p below 1e-15 shown as "p < 1e-15" (no other figure printed a p of 0).
- N-1c Supplementary Fig. S1 legend: "(both mesenteric adipose tissue, GSE123875; …)"; bin filter stated.
- N-1d Cover letter and response letter: "(Fig. 1A) was not affected; the schematic in Fig. 1B, which quoted the F1 values, has been updated".
- N-2a Results and Note 2: "the other eleven exonic SNPs lie outside the part that the pseudogenes copy" (10 lie in exon segments absent from ENSMUST00000146867, the transcript the pseudogenes copy; 17:24,940,701 lies in that transcript just past the copied span).
- N-2b Fig. 3A key: "(n = … tested in liver)" and "catalog genes tested in liver: 117 of 138".
- N-2c TIFFs written as RGB on white (no alpha).
- N-3 every figure re-laid out at its print size, all text ≥ 7 pt (sizes measured from the figure at save time):

| figure | width (mm) | height (mm) | smallest text (pt) |
|---|---|---|---|
| Fig1_survey | 170.0 | 124.5 | 7.0 |
| Fig2_f1_ase | 170.0 | 175.3 | 7.0 |
| Fig3_strain_comparison | 170.0 | 167.6 | 7.0 |
| Fig4_predictors | 170.0 | 185.4 | 7.0 |
| Fig5_catalog | 170.0 | 203.2 | 7.0 |
| Fig6_worked_examples | 170.0 | 193.0 | 7.0 |
| SFig1_GSE55792_haplotype_blocks | 170.0 | 83.8 | 7.0 |
| SFig2_multimapping | 170.0 | 210.8 | 7.0 |
| SFig3_bias_decomposition | 85.0 | 109.2 | 7.0 |

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig1_survey.pdf` | `b4d68c0da3ccab75…` | `f71518378c8c2694…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig1_survey.tiff` | `a6f27841c89d7cc8…` | `b144c239491f8251…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig2_f1_ase.pdf` | `1d0e61c745976d6c…` | `9a3a2b71ce68a470…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig2_f1_ase.tiff` | `827aec957be8766a…` | `7309406fe2c22a40…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig3_strain_comparison.pdf` | `17bb48050479a314…` | `3488226da4a04635…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig3_strain_comparison.tiff` | `6e9fbb9116de1bdf…` | `e10bd3a35a4faad0…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig4_predictors.pdf` | `058002d7aed768e1…` | `46b1beb98265bb2b…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig4_predictors.tiff` | `bd7ac5807062d7f4…` | `d94d6fa186917343…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig5_catalog.pdf` | `e52fd0d47b8a62a0…` | `cf6e7a4462f07aa9…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig5_catalog.tiff` | `cc26be17eb4f07fc…` | `153584967699b44d…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig6_worked_examples.pdf` | `f599bcb13bc0a400…` | `fc2d55a60c146f59…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/Fig6_worked_examples.tiff` | `e5c5476b9ccc45e7…` | `81e7119e6e71c1cf…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/SFig1_GSE55792_haplotype_blocks.pdf` | `74df2005c287e224…` | `3374ea196bf338ce…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/SFig1_GSE55792_haplotype_blocks.tiff` | `0eefe4ff1fd5b09c…` | `a807006a5db1aa29…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/SFig2_multimapping.pdf` | `8c53cb289cd55e29…` | `f8ee604b48cc4bfe…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/SFig2_multimapping.tiff` | `24c8ee10fc3238c3…` | `c1211aefbc9ba092…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/SFig3_bias_decomposition.pdf` | `6b10132da5147391…` | `5ad18b229ed7c588…` |
| N-1a/b, N-2b/c, N-3 | figure re-laid out at print size (RGB TIFF) | `analysis_v2/figures/SFig3_bias_decomposition.tiff` | `77129db5c27dd173…` | `c0256957bf34f25a…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig1_survey.pdf` | `b4d68c0da3ccab75…` | `f71518378c8c2694…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig1_survey.tiff` | `a6f27841c89d7cc8…` | `b144c239491f8251…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig2_f1_ase.pdf` | `1d0e61c745976d6c…` | `9a3a2b71ce68a470…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig2_f1_ase.tiff` | `827aec957be8766a…` | `7309406fe2c22a40…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig3_strain_comparison.pdf` | `17bb48050479a314…` | `3488226da4a04635…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig3_strain_comparison.tiff` | `6e9fbb9116de1bdf…` | `e10bd3a35a4faad0…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig4_predictors.pdf` | `058002d7aed768e1…` | `46b1beb98265bb2b…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig4_predictors.tiff` | `bd7ac5807062d7f4…` | `d94d6fa186917343…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig5_catalog.pdf` | `e52fd0d47b8a62a0…` | `cf6e7a4462f07aa9…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig5_catalog.tiff` | `cc26be17eb4f07fc…` | `153584967699b44d…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig6_worked_examples.pdf` | `f599bcb13bc0a400…` | `fc2d55a60c146f59…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/Fig6_worked_examples.tiff` | `e5c5476b9ccc45e7…` | `81e7119e6e71c1cf…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/SFig1_GSE55792_haplotype_blocks.pdf` | `74df2005c287e224…` | `3374ea196bf338ce…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/SFig1_GSE55792_haplotype_blocks.tiff` | `0eefe4ff1fd5b09c…` | `a807006a5db1aa29…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/SFig2_multimapping.pdf` | `8c53cb289cd55e29…` | `f8ee604b48cc4bfe…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/SFig2_multimapping.tiff` | `24c8ee10fc3238c3…` | `c1211aefbc9ba092…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/SFig3_bias_decomposition.pdf` | `6b10132da5147391…` | `5ad18b229ed7c588…` |
| N-1a/b, N-2b/c, N-3 | figure (submission copy) | `revision_v2/submission/figures/SFig3_bias_decomposition.tiff` | `77129db5c27dd173…` | `c0256957bf34f25a…` |
| N | FIGURES.md5 (18 lines) | `analysis_v2/figures/FIGURES.md5` | `eb46edcbf2ac2c1b…` | `81b7ed8f13cc3428…` |
| N-1d | response letter (Fig. 1A/1B sentence) | `revision_v2/response_to_reviewers.md` | `e03de2a7f4dbecf4…` | `05a03b4d7b559a86…` |
| N-1d | cover letter (Fig. 1A/1B sentence) | `revision_v2/cover_letter.md` | `84979932c300b52f…` | `eaff8de1c070a96a…` |
| N | frozen-input list updated for the figure files | `revision_v2/freeze/inputs.sha256` | `2416d7c113f34bd4…` | `6c045f5b4c935228…` |
| N | deliverables re-sealed | `revision_v2/freeze/deliverables.sha256` | `b3f4dabb7217ce6b…` | `0298c77db005b2f9…` |
| N-1c/d, N-2a | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.docx` | `635d0c3cfc140892…` | `d3cd32b5b31891f5…` |
| N-1c/d, N-2a | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.docx` | `ffb7b6cf9e243783…` | `bf043b6eb9972fa7…` |
| N-1c/d, N-2a | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.docx` | `8aa80783bae22741…` | `9a10e96e5bd76695…` |
| N-1c/d, N-2a | rebuilt deliverable | `revision_v2/submission/cover_letter_v2.docx` | `2d3bfba585b6529b…` | `70230ac56dd43a8f…` |
| N-1c/d, N-2a | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.docx` | `495ca793bfbe652a…` | `7b532f2d52e432e4…` |

## Phase O (2026-10-23): figure label and occlusion fixes — authorised exceptions

Figure-only changes; no data, number or text outside the figures changed: the letters, the xlsx and the .md/.plain.txt sources of both main docx files are byte-identical to Phase N; the two main docx files were re-packed by pandoc during the run (same text, different zip metadata). Every _data.tsv byte-identical; M-1 token differences in the figure PDFs are layout tokens only (axis ticks; Fig. 4D counts now extracted whole); M-2 368 PASS, 0 FAIL.

- O-1a Fig. 3A: only the worked-example genes Ndufs5 and Rps2 labelled, leader lines into the empty upper left, white text background.
- O-1b Fig. 5B: labels stacked beside their column with leader lines, all inside the axes.
- O-1c Supplementary Fig. S2B: the ten labels stacked right of the points with leader lines.
- O-1d Supplementary Fig. S2C: the residual-fraction lists moved right of the bar area (no text on any bar).
- O-1e Supplementary Fig. S2D: horizontal single-letter tissue ticks (A, H, L, S, explained in the legend); one slot per tissue so that the unique and multi-mapping bars of neighbouring tissues no longer overlap.
- O-2a Fig. 2C legend to the upper right (clear of the peak and the mean lines); O-2b Fig. 4A Ndufs5 label in the empty space at x < 0 with a leader line; O-2c Fig. 3B 95th-percentile triangles larger with a black edge; O-2d Fig. 1B solid "maps" arrow from the blue B6-allele box, above the orange box.
- Also found by the automatic check and fixed: horizontal two-line tissue ticks in Figs. 3, 5A and 6 (rotated labels overlapped); Fig. 4D tick labels vertical; Fig. 2 panel letter clear of the tick labels; S2A summary box and legend clear of each other and of the bars.
- O-3 `fig_common.overlap_check` (called by `save()`; a figure is not written if it fails): (i) no two visible text artists overlap (rendered text extents, leader lines excluded); (ii) no text with a background box and no legend frame overlaps a bar or patch of its axes. Exceptions: colour-bar labels and tick labels; texts inside one legend. Every figure was also rendered at 150 dpi and inspected.

| figure | width (mm) | height (mm) | smallest text (pt) | texts checked | overlaps |
|---|---|---|---|---|---|
| Fig1_survey | 170.0 | 124.5 | 7.0 | 37 | 0 |
| Fig2_f1_ase | 170.0 | 175.3 | 7.0 | 89 | 0 |
| Fig3_strain_comparison | 170.0 | 167.6 | 7.0 | 63 | 0 |
| Fig4_predictors | 170.0 | 185.4 | 7.0 | 80 | 0 |
| Fig5_catalog | 170.0 | 203.2 | 7.0 | 198 | 0 |
| Fig6_worked_examples | 170.0 | 193.0 | 7.0 | 88 | 0 |
| SFig1_GSE55792_haplotype_blocks | 170.0 | 83.8 | 7.0 | 34 | 0 |
| SFig2_multimapping | 170.0 | 210.8 | 7.0 | 100 | 0 |
| SFig3_bias_decomposition | 85.0 | 109.2 | 7.0 | 28 | 0 |

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| O | figure labels / occlusion | `analysis_v2/figures/Fig1_survey.pdf` | `2d993ed9c9331c30…` | `f71518378c8c2694…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig1_survey.tiff` | `e0488d254385f0fc…` | `b144c239491f8251…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig2_f1_ase.pdf` | `b6c8dfc513ae7322…` | `9a3a2b71ce68a470…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig2_f1_ase.tiff` | `acf77113d8662289…` | `7309406fe2c22a40…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig3_strain_comparison.pdf` | `f5d38322498edf55…` | `3488226da4a04635…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig3_strain_comparison.tiff` | `a35b5390f302cb52…` | `e10bd3a35a4faad0…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig4_predictors.pdf` | `47cf43576bf3e007…` | `46b1beb98265bb2b…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig4_predictors.tiff` | `081d7cc536acc794…` | `d94d6fa186917343…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig5_catalog.pdf` | `260db644d84b5aad…` | `cf6e7a4462f07aa9…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig5_catalog.tiff` | `8275e0da077959d3…` | `153584967699b44d…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig6_worked_examples.pdf` | `69e6f3b29849d6ad…` | `fc2d55a60c146f59…` |
| O | figure labels / occlusion | `analysis_v2/figures/Fig6_worked_examples.tiff` | `399eff9897b6fd81…` | `81e7119e6e71c1cf…` |
| O | figure labels / occlusion | `analysis_v2/figures/SFig1_GSE55792_haplotype_blocks.pdf` | `5555cd1bb8ec0305…` | `3374ea196bf338ce…` |
| O | figure labels / occlusion | `analysis_v2/figures/SFig2_multimapping.pdf` | `aef1559d0358a236…` | `f8ee604b48cc4bfe…` |
| O | figure labels / occlusion | `analysis_v2/figures/SFig2_multimapping.tiff` | `218fd2b09b35f918…` | `c1211aefbc9ba092…` |
| O | figure labels / occlusion | `analysis_v2/figures/SFig3_bias_decomposition.pdf` | `cc45b9f5340bed88…` | `5ad18b229ed7c588…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig1_survey.pdf` | `2d993ed9c9331c30…` | `f71518378c8c2694…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig1_survey.tiff` | `e0488d254385f0fc…` | `b144c239491f8251…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig2_f1_ase.pdf` | `b6c8dfc513ae7322…` | `9a3a2b71ce68a470…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig2_f1_ase.tiff` | `acf77113d8662289…` | `7309406fe2c22a40…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig3_strain_comparison.pdf` | `f5d38322498edf55…` | `3488226da4a04635…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig3_strain_comparison.tiff` | `a35b5390f302cb52…` | `e10bd3a35a4faad0…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig4_predictors.pdf` | `47cf43576bf3e007…` | `46b1beb98265bb2b…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig4_predictors.tiff` | `081d7cc536acc794…` | `d94d6fa186917343…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig5_catalog.pdf` | `260db644d84b5aad…` | `cf6e7a4462f07aa9…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig5_catalog.tiff` | `8275e0da077959d3…` | `153584967699b44d…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig6_worked_examples.pdf` | `69e6f3b29849d6ad…` | `fc2d55a60c146f59…` |
| O | figure (submission copy) | `revision_v2/submission/figures/Fig6_worked_examples.tiff` | `399eff9897b6fd81…` | `81e7119e6e71c1cf…` |
| O | figure (submission copy) | `revision_v2/submission/figures/SFig1_GSE55792_haplotype_blocks.pdf` | `5555cd1bb8ec0305…` | `3374ea196bf338ce…` |
| O | figure (submission copy) | `revision_v2/submission/figures/SFig2_multimapping.pdf` | `aef1559d0358a236…` | `f8ee604b48cc4bfe…` |
| O | figure (submission copy) | `revision_v2/submission/figures/SFig2_multimapping.tiff` | `218fd2b09b35f918…` | `c1211aefbc9ba092…` |
| O | figure (submission copy) | `revision_v2/submission/figures/SFig3_bias_decomposition.pdf` | `cc45b9f5340bed88…` | `5ad18b229ed7c588…` |
| O | FIGURES.md5 (16 lines) | `analysis_v2/figures/FIGURES.md5` | `5f19e67ec77c17dd…` | `81b7ed8f13cc3428…` |
| O | frozen-input list updated for the figure files | `revision_v2/freeze/inputs.sha256` | `c77b76b556efaa5d…` | `6c045f5b4c935228…` |
| O | deliverables re-sealed | `revision_v2/freeze/deliverables.sha256` | `c965d45012279410…` | `0298c77db005b2f9…` |
| O | docx re-packed by pandoc (text identical: .md and .plain.txt byte-identical) | `revision_v2/submission/LAR_additional_file_1.docx` | `9ab0d10c40b0e638…` | `bf043b6eb9972fa7…` |
| O | docx re-packed by pandoc (text identical: .md and .plain.txt byte-identical) | `revision_v2/submission/LAR_main_text_v2.docx` | `10016ec3e3db77a0…` | `d3cd32b5b31891f5…` |

## Phase P (2026-10-23): answering the reviewers more directly — authorised exceptions

New numbers only through numbers.tsv and M-2 (group P): share of tested genes without an exonic FVB/NJ SNP (48.1%; exons = union of all Ensembl 110 exons of the gene, SNPs = every FVB_NJ SNP record); catalog clusters by a fixed rule (single linkage of catalog genes on one chromosome with Ensembl gene starts < 1 Mb apart, clusters of ≥ 3): 39 genes in 10 clusters, 26 of them in the SNP-induced multi-mapping class; MHC 9 genes in two clusters, replication-dependent histones 11 in two clusters, Mup 3, Ifi200 3. Literature values (30,048 structural variants ≥ 50 bp; the Pde6b nonsense mutation and the 2-bp Hc deletion; ≈50× coverage on mm9) checked against the verbatim sentences of Wong et al. 2012 (revision_v2/audit/wong2012_excerpts.txt). FVB_NJ_v3 checked with NCBI Datasets: GCA_921998635.2, PacBio, falcon-kit 1.8.1, chromosome level, EMBL-EBI, released 2022-07-01, BioProject PRJEB47108 (no publication linked in NCBI; the BioProject is the data deposit of Helmy et al. 2026, Cell Genom 6:101074, whose assemblies are named {strain}_v3). Supplementary Table S2 gained the column cluster; every other cell of every sheet is identical (cell-by-cell comparison). Letters: every quotation of the revised text is checked verbatim against the main text (final_check). P-6 (FVB_NJ_v3 quantification) not run (no approval).

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| P-1b/c | numbers.tsv (new keys only) | `manuscript/revised_v2/generated/numbers.tsv` | `528e8f4b3dddd84e…` | `75b9898203802955…` |
| P-1c | Supplementary Table S2: new column cluster (all other cells identical) | `revision_v2/submission/supplementary_tables_v2.xlsx` | `014b2144c051ef55…` | `b487a26a5933009d…` |
| P-1c | S2 source table | `manuscript/revised_v2/supplementary_tables/TableS2_catalog.tsv` | `873c668665e403ee…` | `82b53a13985bc14a…` |
| P-1b | added: exonic SNP count per gene | `analysis_v2/results/predictors/gene_exonic_snp_counts.tsv` | `(absent)…` | `47e48efc471cc7e8…` |
| P-1c | added: catalog clusters | `analysis_v2/results/catalog/catalog_clusters.tsv` | `(absent)…` | `7c43a063c14b5790…` |
| P-1 | RESULTS.md5 (2 lines added) | `analysis_v2/results/RESULTS.md5` | `f81203f92a161c4c…` | `128502d394c1f9f1…` |
| P-1a | reference list (Helmy 2026 added) | `manuscript/revised_v2/references_source.md` | `a13e27e04f936f20…` | `dbe4270093b207d9…` |
| P-2–P-5 | response letter source | `revision_v2/response_to_reviewers.md` | `aab82cb56f4febf8…` | `05a03b4d7b559a86…` |
| P | frozen-input list | `revision_v2/freeze/inputs.sha256` | `7233afc94d52fc1b…` | `6c045f5b4c935228…` |
| P | deliverables re-sealed | `revision_v2/freeze/deliverables.sha256` | `ca0c33d458e83376…` | `0298c77db005b2f9…` |
| P | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.docx` | `42c4e75482de48bd…` | `bf043b6eb9972fa7…` |
| P | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.md` | `7c610e0a886e9af1…` | `25d1c0c9f05643b8…` |
| P | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.plain.txt` | `557000a03318ac07…` | `70cd3539a84ec502…` |
| P | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.docx` | `7edbb3e47fbd2652…` | `d3cd32b5b31891f5…` |
| P | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.md` | `7d93bb0130f29a3a…` | `7df41fe019e1f212…` |
| P | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.plain.txt` | `2dbc09fd813eca0e…` | `52828cdfa7a7dd35…` |
| P | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.docx` | `f0ba709ade5d5f6f…` | `7b532f2d52e432e4…` |
| P | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.plain.txt` | `0ef423e923c46f2d…` | `4473f0d0e1bda58e…` |
| P | rebuilt deliverable | `revision_v2/submission/cover_letter_v2.docx` | `0de42121ef2aa94c…` | `70230ac56dd43a8f…` |
| P | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.docx` | `5d5a4a09677db35f…` | `9a10e96e5bd76695…` |
| P | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.plain.txt` | `871f18ec443667ce…` | `fbc7827c2bd4002f…` |
| P | rebuilt deliverable | `revision_v2/submission/supplementary_tables_v2.xlsx` | `014b2144c051ef55…` | `b487a26a5933009d…` |

## Phase Q (2026-10-23): claims aligned with the decomposition; FVB_NJ_v3 rationale; within-strain wording; Table 1 — authorised exceptions

Text, letters and the Table 1 layout only; numbers.tsv, the xlsx and every figure unchanged. New text numbers (57 of 138; Ndufs5 bias_matched 1.47–4.61 and b6_loss +1.85 to +5.63; 598 and 539 DE genes; the v1 F1 values in the cover letter) are M-2 claims of group Q. Ensembl annotation of FVB_NJ_v3 verified online (Ensembl release 115, Mus_musculus_fvbnj.FVB_NJ_v3.115.gtf.gz; genebuild last updated 2023-06; own gene identifiers ENSMUSG00205…); cited as reference 'Ensembl FVB 2025'. Table 1: pandoc's cell paragraph style 'Compact' is absent from the localised reference docx, which made LibreOffice collapse the table; cell paragraphs now use the default style, every table has explicit dxa widths, the Table Grid style, 8-pt text and a repeated header row, and Table 1 is in its own landscape section; final_check converts the docx with LibreOffice and checks the layout (freeze/table1_libreoffice_check.png).

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| Q-1c | reference list (Ensembl FVB_NJ_v3 annotation added) | \`manuscript/revised_v2/references_source.md\` | \`2ab63279295afd24…\` | \`dbe4270093b207d9…\` |
| Q-1b/Q-3/Q-4 | response letter source | \`revision_v2/response_to_reviewers.md\` | \`f4e41f6a789c2275…\` | \`05a03b4d7b559a86…\` |
| Q-3/Q-4 | cover letter source | \`revision_v2/cover_letter.md\` | \`62ad198e17fcb464…\` | \`eaff8de1c070a96a…\` |
| Q-5 | assembly: tables with explicit widths, Table 1 landscape | \`revision_v2/scripts/assemble_lar.py\` | \`b4608f140d4ee07e…\` | \`adde9e7901e166b7…\` |
| Q | deliverables re-sealed | \`revision_v2/freeze/deliverables.sha256\` | \`85de2bbf09812085…\` | \`0298c77db005b2f9…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/LAR_additional_file_1.docx\` | \`201f36053e2d9ccd…\` | \`bf043b6eb9972fa7…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/LAR_main_text_v2.docx\` | \`22c349583700fe62…\` | \`d3cd32b5b31891f5…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/LAR_main_text_v2.md\` | \`b4f516b23b07ff18…\` | \`7df41fe019e1f212…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/LAR_main_text_v2.plain.txt\` | \`1f8447950b55bd35…\` | \`52828cdfa7a7dd35…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/LAR_main_text_v2_marked.docx\` | \`31a5c982dd0b5585…\` | \`7b532f2d52e432e4…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/LAR_main_text_v2_marked.plain.txt\` | \`db304961962bce6c…\` | \`4473f0d0e1bda58e…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/cover_letter_v2.docx\` | \`02c89198f62e27eb…\` | \`70230ac56dd43a8f…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/cover_letter_v2.plain.txt\` | \`1d96cad87810359e…\` | \`5d2b2f90bd7daef9…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/response_to_reviewers_v2.docx\` | \`757666198c0a43e2…\` | \`9a10e96e5bd76695…\` |
| Q | rebuilt deliverable | \`revision_v2/submission/response_to_reviewers_v2.plain.txt\` | \`a20a107ea66de3cc…\` | \`fbc7827c2bd4002f…\` |

## Phase R (2026-10-23): three wording corrections and one letter deletion — authorised exceptions

Text only; numbers.tsv, the xlsx and every figure unchanged. 'Pure-strain animals' restricted to the FVB animals of the pure-strain comparison; the 83 animals of the three-arm analysis described as 36 FVB/NJ, 34 C57BL/6J and 13 B6-congenic PyVT tumours (manifest keys; S1b confirms the PyVT tumours were evaluated at the 10% site subset); the negative control described as on a predominantly C57BL/6J background (residual chromosome 16 donor segment); 'detected with similar power' replaced by the overlap of the DE genes; the response-letter sentence referring to an earlier draft deleted. M-2 needle J4 updated to the new wording (value unchanged). cover_letter_v2.docx re-packed by md_to_docx with unchanged source and plain text.

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| R | source script | `revision_v2/scripts/build_phase_j.py` | `7bc9b5b0cf85443d…` | `ba9fbbd6e0741707…` |
| R | source script | `revision_v2/scripts/render_methods_datasets.py` | `5dde78cea6c72f2f…` | `040515f889ca1f76…` |
| R-2 | generated Methods datasets block | `manuscript/revised_v2/generated/methods_datasets.md` | `5827966c75528452…` | `3cad5c8cdadda339…` |
| R-1/R-2/R-4 | response letter source | `revision_v2/response_to_reviewers.md` | `02af69e51f47430a…` | `05a03b4d7b559a86…` |
| R | deliverables re-sealed | `revision_v2/freeze/deliverables.sha256` | `96482d56fedf966e…` | `0298c77db005b2f9…` |
| R | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.docx` | `1a33071ed0800b35…` | `bf043b6eb9972fa7…` |
| R | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.docx` | `1c1462b57d8bc6c6…` | `d3cd32b5b31891f5…` |
| R | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.md` | `308ab6b5706a7dc9…` | `7df41fe019e1f212…` |
| R | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.plain.txt` | `943ff9b1f2d41582…` | `52828cdfa7a7dd35…` |
| R | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.docx` | `bf8ecef60a68d02a…` | `7b532f2d52e432e4…` |
| R | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.plain.txt` | `87ad693ff4854a1e…` | `4473f0d0e1bda58e…` |
| R | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.docx` | `3936481f5858295f…` | `9a10e96e5bd76695…` |
| R | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.plain.txt` | `fd79e5b35beb4a9e…` | `fbc7827c2bd4002f…` |
| R | re-packed by md_to_docx (source and plain text unchanged) | `revision_v2/submission/cover_letter_v2.docx` | `c9c1e6f8afb58e58…` | `70230ac56dd43a8f…` |
| R | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.md` | `0d5a369c84abee55…` | `25d1c0c9f05643b8…` |
| R | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.plain.txt` | `053032cc593efd30…` | `70cd3539a84ec502…` |
| R | checklist regenerated | `revision_v2/submission/submission_checklist.md` | `f9f4d9a1a9a6c021…` | `43e7242ad72753b5…` |

## Phase R2 (2026-10-23): site subset of the identity check — authorised exceptions

Text only (Methods identity sentence; S1b caption site_set); numbers.tsv, the xlsx and every figure unchanged. S1b: the four rejected GSE55788 animals carry site_set sub10pct (20,351 sites); Supplementary Fig. S1 uses their full-site per-site tables. Only files whose bytes changed are listed.

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| R2 | source script | `revision_v2/scripts/build_phase_j.py` | `94cec69e742a570e…` | `ba9fbbd6e0741707…` |
| R2 | deliverables re-sealed | `revision_v2/freeze/deliverables.sha256` | `365f1d96d8b571b3…` | `0298c77db005b2f9…` |
| R2 | rebuilt deliverable | `revision_v2/submission/cover_letter_v2.docx` | `52e72e4c4b3e4d8a…` | `70230ac56dd43a8f…` |
| R2 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.docx` | `eecb64a7e716053c…` | `bf043b6eb9972fa7…` |
| R2 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.md` | `f5fd07ac95d52aa7…` | `25d1c0c9f05643b8…` |
| R2 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.plain.txt` | `39fb256afc48b38c…` | `70cd3539a84ec502…` |
| R2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.docx` | `e2ae1717cd7a8503…` | `d3cd32b5b31891f5…` |
| R2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.docx` | `272959488355b676…` | `7b532f2d52e432e4…` |
| R2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.plain.txt` | `be43d6693c3b1a07…` | `4473f0d0e1bda58e…` |
| R2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.md` | `04c59cd8e7ac0eae…` | `7df41fe019e1f212…` |
| R2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.plain.txt` | `6fa019c283f81abd…` | `52828cdfa7a7dd35…` |
| R2 | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.docx` | `676c49c517a99e5b…` | `9a10e96e5bd76695…` |
| R2 | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.plain.txt` | `0dd060c634b0f990…` | `fbc7827c2bd4002f…` |
| R2 | rebuilt deliverable | `revision_v2/submission/submission_checklist.md` | `d22c93abc51a61bb…` | `43e7242ad72753b5…` |

## Phase S (2026-10-23): ASCII command text in the built docx; survey category definitions — authorised exceptions

Text and build only; numbers.tsv, the xlsx and every figure unchanged. S-1: pandoc smart had turned command options and command quotes into typographic characters in the built docx; assemble_lar.py now escapes them (code spans untouched, prose typography unchanged); final_check checks the built text of all five deliverables (hits before: main 24, marked 24, Additional file 1 2, letters 0; after: 0). S-2: survey categories defined from scripts/stage3_survey_v2.py (STRAIN_AWARE_KW); survey table 140/14/0 of 154, all with sample-level data-processing text. Only deliverables whose bytes changed are listed.

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| S | source | `revision_v2/scripts/build_phase_j.py` | `2400153814c1d60b…` | `ba9fbbd6e0741707…` |
| S | source | `revision_v2/scripts/final_check.sh` | `5ebb60401d821abb…` | `ff75396cec62cc91…` |
| S | source | `revision_v2/response_to_reviewers.md` | `6c38889e66e1b11e…` | `05a03b4d7b559a86…` |
| S-1 | assembly: ASCII command text under pandoc smart | `revision_v2/scripts/assemble_lar.py` | `6e36f9211e267f92…` | `adde9e7901e166b7…` |
| S | deliverables re-sealed | `revision_v2/freeze/deliverables.sha256` | `deb2533d5721ba3f…` | `0298c77db005b2f9…` |
| S | rebuilt deliverable | `revision_v2/submission/cover_letter_v2.docx` | `50e0e4075f98d992…` | `70230ac56dd43a8f…` |
| S | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.docx` | `2b43ede472cdda5b…` | `bf043b6eb9972fa7…` |
| S | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.md` | `ff551e9878c9e72f…` | `25d1c0c9f05643b8…` |
| S | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.plain.txt` | `15539dc63de85ee8…` | `70cd3539a84ec502…` |
| S | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.docx` | `6b5e4ce33643de2c…` | `d3cd32b5b31891f5…` |
| S | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.docx` | `5437a0be491fae92…` | `7b532f2d52e432e4…` |
| S | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.plain.txt` | `6fe4f7ef500ebf15…` | `4473f0d0e1bda58e…` |
| S | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.md` | `9ee3510d5366a6a7…` | `7df41fe019e1f212…` |
| S | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.plain.txt` | `7a6abb8956bb554c…` | `52828cdfa7a7dd35…` |
| S | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.docx` | `60a89e507cc04288…` | `9a10e96e5bd76695…` |
| S | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.plain.txt` | `0dd060c634b0f990…` | `fbc7827c2bd4002f…` |
| S | rebuilt deliverable | `revision_v2/submission/submission_checklist.md` | `017dca62ac709730…` | `43e7242ad72753b5…` |

## Phase S2 (2026-10-23): survey reproducible from the public repository — authorised exceptions

Text and repository only; numbers.tsv, the xlsx and every figure unchanged. survey/ in the repository: survey script (absolute local paths replaced by arguments; no API key or e-mail present), classification table (154; 140/14/0; identical to Supplementary Table S5 in gse/category/keyword), README (full June 2026 GEO texts not cached). token_audit: the 600-character truncation length classified STRUCT. Only deliverables whose bytes changed are listed.

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| S2 | source | `revision_v2/scripts/build_phase_j.py` | `410f84931781b4ba…` | `ba9fbbd6e0741707…` |
| S2 | source | `revision_v2/scripts/commit_phase_h.sh` | `9b5683b07f7f1d01…` | `3b34e95d9684def0…` |
| S2 | source | `revision_v2/scripts/final_check.sh` | `26389460dabbd2bc…` | `ff75396cec62cc91…` |
| S2 | deliverables re-sealed | `revision_v2/freeze/deliverables.sha256` | `d4e65477252b1cd1…` | `0298c77db005b2f9…` |
| S2 | rebuilt deliverable | `revision_v2/submission/cover_letter_v2.docx` | `f2a11bd3751c53f9…` | `70230ac56dd43a8f…` |
| S2 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.docx` | `d8e8e2c0618023bc…` | `bf043b6eb9972fa7…` |
| S2 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.md` | `d56b10871077e6d0…` | `25d1c0c9f05643b8…` |
| S2 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.plain.txt` | `b6aad86c76d47ea8…` | `70cd3539a84ec502…` |
| S2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.docx` | `4cb70e42326e427b…` | `d3cd32b5b31891f5…` |
| S2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.docx` | `b9c960b617319f0b…` | `7b532f2d52e432e4…` |
| S2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.plain.txt` | `f8d00b9f65029f2f…` | `4473f0d0e1bda58e…` |
| S2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.md` | `41547cb9968aabed…` | `7df41fe019e1f212…` |
| S2 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.plain.txt` | `d85ec332590abbc3…` | `52828cdfa7a7dd35…` |
| S2 | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.docx` | `a3b2c737912be072…` | `9a10e96e5bd76695…` |
| S2 | rebuilt deliverable | `revision_v2/submission/submission_checklist.md` | `579fb9594cfc18f1…` | `43e7242ad72753b5…` |

## Phase S3 (2026-10-23): public repository in line with the revised manuscript — authorised exceptions

No analysis, no number changes. Repository README/CITATION generated (numbers from numbers.tsv, in the M-1 audit as 'repo README.md'); v1 scripts moved to legacy/v1/ on the branch; Methods/S5 survey-field sentences corrected to what stage3_survey_v2.py does. Repository files are hashed at HEAD cc9ae4d (before) and in the working tree (after). Only deliverables whose bytes changed are listed.

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| S3 | source | `revision_v2/scripts/build_phase_j.py` | `4b3c4d81f94ba38e…` | `ba9fbbd6e0741707…` |
| S3 | source | `revision_v2/scripts/commit_phase_h.sh` | `501d0b36e8d82005…` | `3b34e95d9684def0…` |
| S3 | source | `revision_v2/scripts/final_check.sh` | `daf7fb8cc610f333…` | `ff75396cec62cc91…` |
| S3 | source | `revision_v2/freeze/token_audit.py` | `ec0d7c5616e219db…` | `fc436f2e707015fe…` |
| S3-1/2 | repository file (HEAD cc9ae4d) | `github_repo/README.md` | `3ee871b29bc1a1f8…` | `9b05e0e131c9301b…` |
| S3-1/2 | repository file (HEAD cc9ae4d) | `github_repo/CITATION.cff` | `b20c6487b1f03bc6…` | `35b467958cd98d93…` |
| S3-1/2 | repository file (HEAD cc9ae4d) | `github_repo/scripts/pipeline/01_build_references.sh` | `3c4cac7371c36275…` | `de0e5e9362806972…` |
| S3 | deliverables re-sealed | `revision_v2/freeze/deliverables.sha256` | `7c4eb64e0cd6e647…` | `0298c77db005b2f9…` |
| S3 | rebuilt deliverable | `revision_v2/submission/cover_letter_v2.docx` | `571d380e01e950a9…` | `70230ac56dd43a8f…` |
| S3 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.docx` | `40fbdcee596bd1e9…` | `bf043b6eb9972fa7…` |
| S3 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.md` | `1550b43383a6b5f4…` | `25d1c0c9f05643b8…` |
| S3 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.plain.txt` | `a5acbfdda1c0fba8…` | `70cd3539a84ec502…` |
| S3 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.docx` | `c54bcc4efc59b275…` | `d3cd32b5b31891f5…` |
| S3 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.docx` | `1e90c6bbde2e6592…` | `7b532f2d52e432e4…` |
| S3 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.plain.txt` | `30b8b516f60e5775…` | `4473f0d0e1bda58e…` |
| S3 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.md` | `ef9fb45219aabc91…` | `7df41fe019e1f212…` |
| S3 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.plain.txt` | `0235e174c67e0282…` | `52828cdfa7a7dd35…` |
| S3 | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.docx` | `b683dda4606dab04…` | `9a10e96e5bd76695…` |
| S3 | rebuilt deliverable | `revision_v2/submission/submission_checklist.md` | `4e05759732334d57…` | `43e7242ad72753b5…` |

## Phase S4 (2026-10-23): clean working tree, identity-check subset, wording — authorised exceptions

No analysis, no number changes. v1-era working-tree copies moved to revision_v2/backup_v1_worktree/; committed modes restored; identity_check/ added (subset reproduced exactly from the documented seed); Methods survey sentence shortened. Only deliverables whose bytes changed are listed.

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| S4 | source | `revision_v2/scripts/build_phase_j.py` | `15bc347175735cc5…` | `ba9fbbd6e0741707…` |
| S4 | source | `revision_v2/scripts/commit_phase_h.sh` | `7aab06700cf249a6…` | `3b34e95d9684def0…` |
| S4 | source | `revision_v2/scripts/final_check.sh` | `1f0954cc7a4a0fe4…` | `ff75396cec62cc91…` |
| S4 | source | `revision_v2/scripts/make_repo_readme.py` | `200fa9c9c9a34527…` | `812a269f56f5ebba…` |
| S4 | repository file (HEAD fbff905) | `github_repo/README.md` | `a052fa5ed64940c3…` | `9b05e0e131c9301b…` |
| S4 | deliverables re-sealed | `revision_v2/freeze/deliverables.sha256` | `a6d472f356bc5e65…` | `0298c77db005b2f9…` |
| S4 | rebuilt deliverable | `revision_v2/submission/cover_letter_v2.docx` | `9397051f429b1709…` | `70230ac56dd43a8f…` |
| S4 | rebuilt deliverable | `revision_v2/submission/LAR_additional_file_1.docx` | `08e8ce52035cde7f…` | `bf043b6eb9972fa7…` |
| S4 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.docx` | `0f230370284b26da…` | `d3cd32b5b31891f5…` |
| S4 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.docx` | `594fee1fc2e31006…` | `7b532f2d52e432e4…` |
| S4 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2_marked.plain.txt` | `305448b8d2a67b75…` | `4473f0d0e1bda58e…` |
| S4 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.md` | `885bf82236b9c3f4…` | `7df41fe019e1f212…` |
| S4 | rebuilt deliverable | `revision_v2/submission/LAR_main_text_v2.plain.txt` | `5f0065cfb15f0424…` | `52828cdfa7a7dd35…` |
| S4 | rebuilt deliverable | `revision_v2/submission/response_to_reviewers_v2.docx` | `b6d0a41694d53365…` | `9a10e96e5bd76695…` |
| S4 | rebuilt deliverable | `revision_v2/submission/submission_checklist.md` | `738dc818e4870005…` | `43e7242ad72753b5…` |

## Phase S5 (2026-10-24): CITATION.cff only — authorised exceptions

No analysis, no number changes, manuscript text unchanged. CITATION.cff (generated by make_repo_readme.py): preferred-citation status: submitted; author e-mail removed; date-released 2026-10-24 (planned push date); valid against CFF 1.2.0 (cffconvert 2.0.0 --validate). Sources are hashed at HEAD 0c967eb (before) and in the working tree (after).

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| S5 | source | `revision_v2/scripts/make_repo_readme.py` | `e5e8b14873be226d…` | `812a269f56f5ebba…` |
| S5 | source | `revision_v2/scripts/final_check.sh` | `e0529fe973e2ae97…` | `ff75396cec62cc91…` |
| S5 | source | `revision_v2/freeze/make_freeze_md.py` | `02467968471534b0…` | `572f9f2f8d31c7f1…` |
| S5 | source | `revision_v2/CHANGELOG.md` | `3938a65940522e2c…` | `169e81384cb75f28…` |
| S5 | repository file | `repo_release/CITATION.cff` | `35b467958cd98d93…` | `603f34d130acea0b…` |
