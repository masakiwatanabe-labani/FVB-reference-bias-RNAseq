# Phase F — final check (2026-10-20 09:47 JST)

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
- [x] Appendix C grep: PASS
- (exempt, data tables) EXEMPT (data) [cardiac]: 2 line(s): manuscript/revised_v2/supplementary_tables/TableS5_survey_of_practice.tsv 
- (exempt, data tables) EXEMPT (data) [heart]: 3 line(s): manuscript/revised_v2/supplementary_tables/TableS5_survey_of_practice.tsv 
- [x] every [CITE] key (surname + year) has a reference entry
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
- [x] manuscript/revised_v2/discussion.md (2652 words)
- [x] manuscript/revised_v2/methods.md (3904 words)
- [x] manuscript/revised_v2/figure_legends.md (1835 words)
- [x] manuscript/revised_v2/table1.md (410 words)
- [x] manuscript/revised_v2/supplementary_table_captions.md (1114 words)
- [x] manuscript/revised_v2/supplementary_note1.md (797 words)
- [x] manuscript/revised_v2/supplementary_note2.md (1115 words)
- [x] manuscript/revised_v2/references.md (634 words)
- [x] manuscript/revised_v2/manuscript_combined.md (17329 words)
- [x] revision_v2/text_changes.md (1768 words)
- [x] revision_v2/response_to_reviewers.md (2232 words)
- [x] revision_v2/cover_letter.md (512 words)
- [x] revision_v2/summary_of_changes.md (452 words)
- [x] revision_v2/CHANGELOG.md (8487 words)

## 6. Open items for the authors (human decision)
- Appendix C, decided 2026-10-20: '+1.14' (Rps2 liver bias, numbers.tsv ex_extreme_liver_bias) and '0.005' (liver WASP − naive mean, gw_liver_wasp_minus_naive_mean_logFC) are coincidental matches of NEW values with v1 strings; new source verified in numbers.tsv, values kept as they are. 'heart' is the GEO tissue attribute of the rejected GSE55792 samples; prose and tables unified to 'heart'. These contexts are ACCEPTED in appendix_c_grep.sh.
- Doran citation, decided 2026-10-20: [CITE Doran 2021] replaced by Doran et al. 2016 (Genome Biol 17:167, doi:10.1186/s13059-016-1024-y; DOI resolves — doi.org 302 to Genome Biology, Crossref record retrieved) as the MGP variant-catalogue reference; entry 13 in references.md.
- Discussion states that the FVB_NJ_v3 comparison (optional Phase G) was not performed.
- github_repo has uncommitted v1-era modifications (README.md, CITATION.cff, LICENSE, environment.yml, .gitignore); prepare_github_v2.sh commits only the v2 files it adds plus the README paragraph — review before push.

**Result: 0 failing check(s).**
