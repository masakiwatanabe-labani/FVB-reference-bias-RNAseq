# Number freeze — Phase M (2026-10-23)

Sealed in the commit that adds this file (parent 2335065, branch v2-corrected); that commit is the one the authors tag v2-corrected-r1.

## What was verified

- **M-0 inputs frozen**: `inputs.sha256` lists every frozen input (numbers.tsv, supplementary_tables_v2.xlsx, every file of RESULTS.md5 and FIGURES.md5, and the two md5 lists); final_check fails if any changes. `build_submission.sh` does not rewrite the xlsx or the figures once frozen.
- **M-1 numeric-token audit** (`token_audit.py`, `resolve_tokens.py` → `token_audit.tsv`): 1578 numeric tokens in the main text (all sections, Table 1, legends), Additional file 1, both letters and the text of the nine figure PDFs.
- **M-2 independent recomputation** (`recompute_independent.py` → `recompute_report.tsv`): no import of the pipeline or of the text/number generators, no numbers.tsv; primary result files only (edgeR tables, count matrices and featureCounts summaries, STAR junction lists and Log.out, F1 and identity per-site tables, WASP allele-removal table, within-strain DE tables, FVB_NJ VCFs and the Ensembl 110 GTF via bcftools/bedtools, rps2_* files, vW tag counts, survey table, submitted-version catalog for the letters). Bounds pass only if the exact value satisfies the bound as printed.
- **D3 double-rounding scan** (`double_rounding_scan.py`): every key printed in the text compared with the value rounded once from full precision (`generated/numbers_fullprecision.tsv`).
- **Seal**: `deliverables.sha256` (every file in revision_v2/submission/, figures included).

## Token classes (M-1)

| class | tokens |
|---|---|
| STRUCT | 793 |
| KEY | 773 |
| M2 | 9 |
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
| **total** | **354** | **0** | **2** |

## Authorised exceptions to the freeze (authors' decisions, 2026-10-23)

| id | change | file | sha256 before | sha256 after |
|---|---|---|---|---|
| D1 | Fig. 1 regenerated (v2 F1 table) | `analysis_v2/figures/Fig1_survey.pdf` | `2c37ea19711eda5b…` | `b4d68c0da3ccab75…` |
| D1 | Fig. 1 regenerated | `analysis_v2/figures/Fig1_survey.tiff` | `c10dc16e6fbc38b9…` | `a6f27841c89d7cc8…` |
| D1 | Fig. 1 (submission copy) | `revision_v2/submission/figures/Fig1_survey.pdf` | `2c37ea19711eda5b…` | `b4d68c0da3ccab75…` |
| D1 | Fig. 1 (submission copy) | `revision_v2/submission/figures/Fig1_survey.tiff` | `c10dc16e6fbc38b9…` | `a6f27841c89d7cc8…` |
| D1b | Suppl. Fig. S1 regenerated (v2 reference animals) | `analysis_v2/figures/SFig1_GSE55792_haplotype_blocks.pdf` | `e5d6642cc519cead…` | `74df2005c287e224…` |
| D1b | Suppl. Fig. S1 regenerated | `analysis_v2/figures/SFig1_GSE55792_haplotype_blocks.tiff` | `c68a08777c22c2b6…` | `0eefe4ff1fd5b09c…` |
| D1b | Suppl. Fig. S1 data | `analysis_v2/figures/SFig1_GSE55792_haplotype_blocks_data.tsv` | `8b46b0a4a84629da…` | `6e21a0919d604a90…` |
| D1b | Suppl. Fig. S1 (submission copy) | `revision_v2/submission/figures/SFig1_GSE55792_haplotype_blocks.pdf` | `e5d6642cc519cead…` | `74df2005c287e224…` |
| D1b | Suppl. Fig. S1 (submission copy) | `revision_v2/submission/figures/SFig1_GSE55792_haplotype_blocks.tiff` | `c68a08777c22c2b6…` | `0eefe4ff1fd5b09c…` |
| D1b | added: reference-animal candidates (+6 per-site tables) | `analysis_v2/results/identity/full_site_tables/candidates.tsv` | `(absent)…` | `0a164012e848db18…` |
| D2 | added: Rps2 exonic variants, union definition | `analysis_v2/results/worked_examples/rps2_exonic_variants_union.tsv` | `(absent)…` | `fbf5d346b6f3cbc5…` |
| D2 | unchanged (no new SNP inside the copied part) | `analysis_v2/results/worked_examples/rps2_pseudogene_alleles.tsv` | `8cf1d351f6a169a0…` | `8cf1d351f6a169a0…` |
| D2/D3/D4/D5/E1 | numbers.tsv (3 values, 59 CPM precisions, new keys) | `manuscript/revised_v2/generated/numbers.tsv` | `b7b3af78ef8afd3f…` | `528e8f4b3dddd84e…` |
| D1/D1b | FIGURES.md5 (5 lines) | `analysis_v2/figures/FIGURES.md5` | `a15db3c5e3fb67f1…` | `eb46edcbf2ac2c1b…` |
| D1b/D2 | RESULTS.md5 (8 lines added) | `analysis_v2/results/RESULTS.md5` | `402e3f9f84bfc3d5…` | `f81203f92a161c4c…` |
| D5b | Table 1 tissue cell | `manuscript/revised_v2/table1.md` | `71ab84cee173b941…` | `dd305f1db254fb99…` |
| M-0 | frozen-input list updated for the files above | `revision_v2/freeze/inputs.sha256` | `351616284102d0a4…` | `2416d7c113f34bd4…` |

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
