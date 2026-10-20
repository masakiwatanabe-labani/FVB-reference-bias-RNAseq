# Phase F — final check (2026-10-20 09:40 JST)

## 1. Integrity of results and figures
- [x] analysis_v2/results: all 161 TSV files match RESULTS.md5
- [x] analysis_v2/figures: all 24 files match FIGURES.md5 (Fig1–6, SFig1, SFig2: PDF + TIFF + _data.tsv)
- [x] 83 done.all markers (83 animals processed through the three arms)
- [x] no STOP flag; scheduler log ends with COMPLETE: 2 line(s)

## 2. Manifest-driven generation (rule 0.4-1): re-render and diff
- [x] render_methods_datasets.py ran
- [x] make_table1.py ran
- [x] table1.md == generated/table1.md
- [x] Methods 'RNA-seq datasets' block == generated/methods_datasets.md (BEGIN/END markers)
- [x] supplementary tables rebuild byte-identical (12 sheets)
- [x] results_numbers.py ran: 593 numbers with evidence paths

## 3. Manuscript checks
- [x] figure citation order: all figures and panels in numerical order (check_figure_order.py --legend-check)
- [ ] **Appendix C grep: FAILED** — hits:
      FORBIDDEN [+1.14]:
      manuscript/revised_v2/manuscript_combined.md:97:*Rps2* (ribosomal protein S2; chromosome 17; 53 FVB/NJ SNPs in a 3.8-kb gene body, 13.9 SNPs/kb) illustrates the
      manuscript/revised_v2/results.md:67:*Rps2* (ribosomal protein S2; chromosome 17; 53 FVB/NJ SNPs in a 3.8-kb gene body, 13.9 SNPs/kb) illustrates the direct-misa
      FORBIDDEN [0.005]:
      manuscript/revised_v2/manuscript_combined.md:61:Reflecting this minimal read-level attrition, the genome-wide mean log₂FC under WASP alignment differed from t
      manuscript/revised_v2/results.md:31:Reflecting this minimal read-level attrition, the genome-wide mean log₂FC under WASP alignment differed from the naive arm
- (exempt, data tables) EXEMPT (data) [cardiac]: 2 line(s): manuscript/revised_v2/supplementary_tables/TableS5_survey_of_practice.tsv 
- (exempt, data tables) EXEMPT (data) [heart]: 60 line(s): analysis_v2/figures/SFig1_GSE55792_haplotype_blocks_data.tsv manuscript/revised_v2/supplementary_tables/TableS1a_manifest.tsv manuscript/revised_v2/supplementary_tables/TableS1b_identity_verification.tsv manuscript/revised_v2/supplementary_tables/TableS5_survey_of_practice.tsv 
- [ ] **FAIL** [CITE] keys without a reference entry: Doran 2021
- [x] Availability names branch/tag v2-corrected

## 4. Rules 0.4
- [x] 0.4-1 manifest first: Table 1, Methods datasets, counts/sample sheets (merge_counts.py reads manifest.tsv only) — evidence: this section 2
- [x] 0.4-2 no SRR block selection: scheduler.sh sel() and merge_counts.py select by manifest attributes (role, include, tissue, strain, design)
- [x] 0.4-3 evidence: numbers.tsv (key → evidence file); CHANGELOG entries name their files
- [x] 0.4-4 identity QC on all samples: qc/sample_identity_verification.tsv 99 rows, 95 used all match, 4 failed listed (Supplementary Table S1b)
- [x] 0.4-5 anomalies → identity re-check before explanation: GSE55792 (audit/gse55792_qc_note.md), PyVT chr 3/10/16 full-site recount (followup_analysis.py section 1)
- [x] 0.4-6 Appendix C grep: see section 3 (residual items are listed in text_changes.md for human decision)
- [x] 0.4-7 legacy isolation: analysis_v2 code reads nothing under legacy/ — 3 matching files (old_vs_new_overlap.py reads legacy on purpose for the response letter only)
- [x] 0.4-8 criteria fixed: params.yaml; catalog 0.5 log2 / ≥ 2 tissues / FVB GTF; predictor definition reverted to old Methods (CHANGELOG 2026-10-19 22:36)
- [x] 0.4-9 GitHub: prepare_github_v2.sh stages branch v2-corrected locally (no push); README paragraph in revision_v2/README_correction_paragraph.md
- [x] 0.4-10 heavy jobs backgrounded with logs: analysis_v2/logs/*.log

## 5. Deliverables
- [x] manuscript/revised_v2/abstract.md (424 words)
- [x] manuscript/revised_v2/introduction.md (991 words)
- [x] manuscript/revised_v2/results.md (3442 words)
- [x] manuscript/revised_v2/discussion.md (2585 words)
- [x] manuscript/revised_v2/methods.md (3904 words)
- [x] manuscript/revised_v2/figure_legends.md (1835 words)
- [x] manuscript/revised_v2/table1.md (410 words)
- [x] manuscript/revised_v2/supplementary_table_captions.md (1114 words)
- [x] manuscript/revised_v2/supplementary_note1.md (797 words)
- [x] manuscript/revised_v2/supplementary_note2.md (1115 words)
- [x] manuscript/revised_v2/references.md (596 words)
- [x] manuscript/revised_v2/manuscript_combined.md (17224 words)
- [x] revision_v2/text_changes.md (1655 words)
- [x] revision_v2/response_to_reviewers.md (2232 words)
- [x] revision_v2/cover_letter.md (512 words)
- [x] revision_v2/summary_of_changes.md (452 words)
- [x] revision_v2/CHANGELOG.md (7952 words)

## 6. Open items for the authors (human decision)
- Appendix C: '+1.14' (Rps2 liver bias) and '0.005' (liver WASP − naive mean) in Results are new values that coincide with forbidden v1 strings; 'heart' = GEO tissue attribute of the rejected GSE55792 samples in Supplementary Table S1a/b and Supplementary Fig. S1 labels (main text says 'myocardial').
- [CITE Doran 2021] has no entry in references.md (pre-existing in v1; DOI not verified here) — add or remove.
- Discussion states that the FVB_NJ_v3 comparison (optional Phase G) was not performed.
- github_repo has uncommitted v1-era modifications (README.md, CITATION.cff, LICENSE, environment.yml, .gitignore); prepare_github_v2.sh commits only the v2 files it adds plus the README paragraph — review before push.

**Result: 2 failing check(s).**
