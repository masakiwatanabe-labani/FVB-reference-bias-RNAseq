# Phase J — final check (2026-10-21 12:08 JST)

## 1. Integrity of results and figures
- [x] analysis_v2/results: all 233 files match RESULTS.md5 (includes f1_ase/v2/: 24 + 24 per-site tables, per-animal/summary tables, catalog decomposition)
- [x] analysis_v2/figures: all 28 files match FIGURES.md5 (Fig 1–6, SFig 1–3: PDF + TIFF + _data.tsv)
- [x] Fig. 2 (regenerated from f1_ase/v2/, with panel D) and Supplementary Fig. S3 are listed in FIGURES.md5
- [x] 83 done.all markers (83 animals processed through the three arms); 12 done.f1_ase markers: 12
- [x] no STOP flag

## 2. Manifest-driven generation (rule 0.4-1): re-render and diff
- [x] render_methods_datasets.py ran
- [x] make_table1.py ran
- [x] table1.md == generated/table1.md
- [x] Methods 'RNA-seq datasets' block == generated/methods_datasets.md (BEGIN/END markers)
- [x] supplementary tables rebuild byte-identical (12 sheets)
- [x] results_numbers.py ran: 732 numbers with evidence paths

## 3. Manuscript checks
- [x] figure citation order (incl. Supplementary Fig. S3): all figures and panels in numerical order (check_figure_order.py --legend-check)
- [x] Appendix C grep (incl. the H-1d additions 0.333, 0.386, 17 percentage, 274, 872): PASS
- (exempt, data tables) EXEMPT (data) [cardiac]: 2 line(s): manuscript/revised_v2/supplementary_tables/TableS5_survey_of_practice.tsv 
- (exempt, data tables) EXEMPT (data) [heart]: 3 line(s): manuscript/revised_v2/supplementary_tables/TableS5_survey_of_practice.tsv 
- [x] no [CITE placeholder left in any main-text or Additional-file section
- [x] numbered citations: 22 references, numbered 1..22 consecutively; all 22 cited; first appearances in increasing order: True
- [x] all DOIs in references.md resolve (HTTP 301/302 from doi.org; revision_v2/audit/doi_check_final.txt)
- [x] structured Abstract: 278 words (≤ 350), 4 headings (Background/Methods/Results/Conclusions)
- [x] assemble_lar.py rebuilt LAR_main_text_v2.docx and LAR_additional_file_1.docx (pandoc, --reference-doc 'LAR text submit.docx')
- [x] section order in the docx (plain dump): Abstract (structured) → Keywords → Background>Methods>Results>Discussion>Conclusions>List of abbreviations>Declarations>References>Figure legends> Table 1
- [x] H-5 grep (corrected analysis / this revision / earlier version / submitted version / reviewer) in the main text and Additional file 1: 0 hits
- [x] all 58 f1_* numbers come from the new per-site tables (evidence column: analysis_v2/results/f1_ase/v2/ or the H-1a audit files); the v1 tables f1_ase/f1_ase_per_sample.tsv are not an evidence source
- [x] docx → plain: no [CITE, no local paths, no continued numbered lists (0 Step paragraphs; the reconstruction keeps the v1 prose steps), Fig. 2 significance asterisks rendered literally (1 line)
- [x] Availability (Declarations) names branch/tag v2-corrected

## 3b. Phase J checks (minimal-change reconstruction on the submitted text)
- [x] change ledger rebuilt: total {'KEEP': 512, 'NEW': 465, 'EDIT': 151, 'DROP': 311}; body sentences 726: KEEP 387 (53%)
- [x] ledger completeness: 1439 rows; 726 body sentences with a row; KEEP sentences not found verbatim in v1: 0; EDIT/NEW numbers not in numbers.tsv: 0; EDIT/NEW/DROP rows without trigger: 0
- [x] response letter: every 'Changed:' pointer resolves to ledger rows with the matching trigger (18 pointers)
- [x] Abstract ≤ 350 words both ways (body 278, with headings 283)
- [x] '3Rs' is used in the text (Discussion), so the abbreviation entry is used
- [x] abbreviation list: 29 entries, unused: []
- [x] Phase I wording present: direct-class read fate (4 places), 'expected to be heterozygous' (Methods), 'Mechanism 1./2.' labels (Note 2)
- [x] J-0 semantic ban list (also in appendix_c_grep.sh): 0 hits in the docx dumps
- [x] marked copy present (77386 bytes; generated from the ledger by mark_from_ledger.py)

## 4. Rules 0.4
- [x] 0.4-1 manifest first: Table 1, Methods datasets, counts/sample sheets (merge_counts.py reads manifest.tsv only) — evidence: section 2
- [x] 0.4-2 no SRR block selection: scheduler.sh sel() and merge_counts.py select by manifest attributes (role, include, tissue, strain, design)
- [x] 0.4-3 evidence: numbers.tsv (key → evidence file); CHANGELOG entries name their files
- [x] 0.4-4 identity QC on all samples: qc/sample_identity_verification.tsv 99 rows, 95 used all match (F1 rows re-checked on the regenerated BAMs, H-1), 4 failed listed (Supplementary Table S1b)
- [x] 0.4-5 anomalies → identity re-check before explanation: GSE55792 (audit/gse55792_qc_note.md), PyVT chr 3/10/16 full-site recount, F1 site-set artefact (audit/f1_old_site_genotypes.md; old-site reproduction f1_ase_summary_v2.oldsites.tsv)
- [x] 0.4-6 Appendix C grep: section 3
- [x] 0.4-7 legacy isolation: analysis_v2 code reads nothing under legacy/ — 4 matching files (old_vs_new_overlap.py reads legacy on purpose for the response letter only)
- [x] 0.4-8 criteria fixed: params.yaml; catalog 0.5 log2 / ≥ 2 tissues / FVB GTF unchanged (138 genes); decomposition added as a description, not a new criterion
- [x] 0.4-9 GitHub: branch v2-corrected committed locally (commit_phase_h.sh); no push; tag to be created by the authors at submission
- [x] 0.4-10 heavy jobs backgrounded with logs: analysis_v2/logs/*.log (f1_ase_rerun.log, f1_oldsites_recount.log)

## 5. Deliverables (revision_v2/submission/)
- [x] revision_v2/submission/LAR_main_text_v2.docx (75873 bytes)
- [x] revision_v2/submission/LAR_additional_file_1.docx (47148 bytes)
- [x] revision_v2/submission/supplementary_tables_v2.xlsx (158097 bytes)
- [x] revision_v2/submission/response_to_reviewers_v2.docx (43894 bytes)
- [x] revision_v2/submission/cover_letter_v2.docx (38540 bytes)
- [x] figures/: 9 PDF + 9 TIFF (Fig 1–6, SFig S1–S3)
- [x] manuscript/revised_v2/abstract_structured.md (284 words)
- [x] manuscript/revised_v2/introduction.md (769 words)
- [x] manuscript/revised_v2/results.md (4564 words)
- [x] manuscript/revised_v2/discussion.md (3294 words)
- [x] manuscript/revised_v2/methods.md (4387 words)
- [x] manuscript/revised_v2/figure_legends.md (2041 words)
- [x] manuscript/revised_v2/table1.md (410 words)
- [x] manuscript/revised_v2/supplementary_table_captions.md (1325 words)
- [x] manuscript/revised_v2/supplementary_note1.md (821 words)
- [x] manuscript/revised_v2/supplementary_note2.md (1255 words)
- [x] manuscript/revised_v2/references.md (591 words)
- [x] manuscript/revised_v2/declarations.md (259 words)
- [x] manuscript/revised_v2/abbreviations.md (132 words)
- [x] manuscript/revised_v2/keywords.md (15 words)
- [x] revision_v2/text_changes.md (3293 words)
- [x] revision_v2/response_to_reviewers.md (2706 words)
- [x] revision_v2/cover_letter.md (611 words)
- [x] revision_v2/CHANGELOG.md (10752 words)

## 6. Open items for the authors (human decision)
- Appendix C, decided 2026-10-20: '+1.14' (Rps2 liver bias) and '0.005' (liver WASP − naive mean) are coincidental matches of NEW values with v1 strings (ACCEPTED contexts in appendix_c_grep.sh); 'heart' is the GEO tissue attribute of the rejected GSE55792 samples.
- The FVB_NJ_v3 comparison (optional Phase G) was not performed; the Discussion says so and names the Rps2 skin locus as the first place to look.
- github_repo has uncommitted v1-era modifications (.gitignore, CITATION.cff, LICENSE, environment.yml, scripts/figures/*); commit_phase_h.sh commits only the v2 files — review before push.

**Result: 0 failing check(s).**
