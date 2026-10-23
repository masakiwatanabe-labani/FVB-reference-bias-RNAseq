# Phase J — final check (2026-10-23 22:02 JST)

## 1. Integrity of results and figures
- [x] analysis_v2/results: all 245 files match RESULTS.md5 (includes f1_ase/v2/: 24 + 24 per-site tables, per-animal/summary tables, catalog decomposition)
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
- [x] results_numbers.py ran: 817 numbers with evidence paths

## 3. Manuscript checks
- [x] figure citation order (incl. Supplementary Fig. S3): all figures and panels in numerical order (check_figure_order.py --legend-check)
- [x] Appendix C grep (incl. the H-1d additions 0.333, 0.386, 17 percentage, 274, 872): PASS
- (exempt, data tables) EXEMPT (data) [cardiac]: 2 line(s): manuscript/revised_v2/supplementary_tables/TableS5_survey_of_practice.tsv 
- (exempt, data tables) EXEMPT (data) [heart]: 3 line(s): manuscript/revised_v2/supplementary_tables/TableS5_survey_of_practice.tsv 
- [x] no [CITE placeholder left in any main-text or Additional-file section
- [x] numbered citations: 24 references, numbered 1..24 consecutively; all 24 cited; first appearances in increasing order: True
- [x] all DOIs in references.md resolve (HTTP 301/302 from doi.org; revision_v2/audit/doi_check_final.txt)
- [x] structured Abstract: 338 words (≤ 350), 4 headings (Background/Methods/Results/Conclusions)
- [x] deliverables sealed (Phase M-4): assemble_lar.py not re-run, because pandoc would rewrite the docx bytes; the sealed files are checked against deliverables.sha256 below
- [x] section order in the docx (plain dump): Abstract (structured) → Keywords → Background>Methods>Results>Discussion>Conclusions>List of abbreviations>Declarations>References>Figure legends> Table 1
- [x] H-5 grep (corrected analysis / this revision / earlier version / submitted version / reviewer) in the main text and Additional file 1: 0 hits
- [x] all 58 f1_* numbers come from the new per-site tables (evidence column: analysis_v2/results/f1_ase/v2/ or the H-1a audit files); the v1 tables f1_ase/f1_ase_per_sample.tsv are not an evidence source
- [x] docx → plain: no [CITE, no local paths, no continued numbered lists (0 Step paragraphs; the reconstruction keeps the v1 prose steps), Fig. 2 significance asterisks rendered literally (1 line)
- [x] Availability (Methods, Declarations) and response item 5 name branch v2-corrected, tag v2-corrected-r1; 'branch and tag' 0 hits (Phase L-2)

## 3b. Phase J checks (minimal-change reconstruction on the submitted text)
- [x] change ledger rebuilt: total {'KEEP': 479, 'NEW': 526, 'EDIT': 166, 'DROP': 333}; body sentences 746: KEEP 361 (48%)
- [x] ledger completeness: 1504 rows; 746 body sentences with a row; KEEP sentences not found verbatim in v1: 0; EDIT/NEW numbers not in numbers.tsv: 0; EDIT/NEW/DROP rows without trigger: 0
- [x] response letter: every 'Changed:' pointer resolves to ledger rows with the matching trigger (21 pointers)
- [x] Abstract ≤ 350 words both ways (body 338, with headings 343)
- [x] '3Rs' is used in the text (Discussion), so the abbreviation entry is used
- [x] abbreviation list: 31 entries, unused: []
- [x] Phase I/K-1′ wording present: read fate 'multi-mappers on both references' (8 places), classes defined operationally by |Δ| (Methods, Note 2), 'expected to be heterozygous' (Methods), 'Mechanism 1./2.' labels (Note 2)
- [x] K-1′ withdrawn read-fate phrases: 0 hits in main text, Additional file 1 and both letters
- [x] every Rps2-skin recovery statement carries 'under unique-only counting'
- [x] Phase L-1 greps: 0 hits for the six withdrawn indel/SNP-only phrases; pseudogene-allele explanation present in Results and Note 2
- [x] all Rps2/pseudogene CPM keys in numbers.tsv use the single definition (evidence rps2_multimap_cpm.tsv)
- [x] J-0 semantic ban list (also in appendix_c_grep.sh): 0 hits in the docx dumps
- [x] every 'largest single-tissue bias' carries the class qualifier (8 occurrences)
- [x] marked copy present (80247 bytes; generated from the ledger by mark_from_ledger.py)
- [x] Phase M freeze: all 277 frozen inputs (numbers.tsv, supplementary_tables_v2.xlsx, RESULTS.md5 and FIGURES.md5 files) match revision_v2/freeze/inputs.sha256
- [x] Phase M-1 token audit: 1784 numeric tokens; UNMAPPED only 1,047 (search hits, no log) and 290 (GEO series size, verified against the GEO SOFT record)
- [x] Phase M-2 independent recomputation: 399 claims PASS, 0 FAIL (recompute_report.tsv)
- [x] Phase M-3 D3 double-rounding scan: keys printed in the text: 221; checked against full precision: 220; double-rounding defects: 0
- [x] Phase M-4 seal: all 34 files of revision_v2/submission/ match deliverables.sha256
- [x] response letter: every quoted passage of the revised text occurs verbatim in the main text (12 quotes, 0 not found)
- [x] ASCII command text in every built docx (no '–option', no curly quote after '=' or in 'samtools view -e'): LAR_main_text_v2 0; LAR_main_text_v2_marked 0; LAR_additional_file_1 0; response_to_reviewers_v2 0; cover_letter_v2 0
- [x] survey release (survey_release/ → survey/): 154 series, 140/14/0 stated_B6/not_reported/strain_aware, table = Supplementary Table S5: True, local paths/e-mail/API keys in script: 0
- [x] repository README/CITATION: v1 strings 0; README accessions 8 = Table 1 accessions 8: True; committed copies on v2-corrected match: True
- [x] identity_check release: 20351 sites = subset VCF used: True; intervals identical: True; committed copies match: True
- [x] github_repo working tree clean (git status --porcelain: 0 lines)
- [x] Table 1 survives PDF conversion (LibreOffice headless): page 34 of 34 (landscape); 7/7 data rows with every cell starting on the row's first line
- [x] reference 'Ensembl FVB 2025' URL resolves (HTTP 200; Ensembl release 115 FVB_NJ_v3 GTF directory); DOIs: revision_v2/audit/doi_check*.txt

## 4. Rules 0.4
- [x] 0.4-1 manifest first: Table 1, Methods datasets, counts/sample sheets (merge_counts.py reads manifest.tsv only) — evidence: section 2
- [x] 0.4-2 no SRR block selection: scheduler.sh sel() and merge_counts.py select by manifest attributes (role, include, tissue, strain, design)
- [x] 0.4-3 evidence: numbers.tsv (key → evidence file); CHANGELOG entries name their files
- [x] 0.4-4 identity QC on all samples: qc/sample_identity_verification.tsv 99 rows, 95 used all match (F1 rows re-checked on the regenerated BAMs, H-1), 4 failed listed (Supplementary Table S1b)
- [x] 0.4-5 anomalies → identity re-check before explanation: GSE55792 (audit/gse55792_qc_note.md), PyVT chr 3/10/16 full-site recount, F1 site-set artefact (audit/f1_old_site_genotypes.md; old-site reproduction f1_ase_summary_v2.oldsites.tsv)
- [x] 0.4-6 Appendix C grep: section 3
- [x] 0.4-7 legacy isolation: analysis_v2 code reads nothing under legacy/ — 4 matching files (old_vs_new_overlap.py reads legacy on purpose for the response letter only)
- [x] 0.4-8 criteria fixed: params.yaml; catalog 0.5 log2 / ≥ 2 tissues / FVB GTF unchanged (138 genes); decomposition added as a description, not a new criterion
- [x] 0.4-9 GitHub: branch v2-corrected committed locally (commit_phase_h.sh); no push; tag v2-corrected-r1 to be created by the authors on the final commit (not v2-corrected: it would clash with the branch name)
- [x] 0.4-10 heavy jobs backgrounded with logs: analysis_v2/logs/*.log (f1_ase_rerun.log, f1_oldsites_recount.log)

## 5. Deliverables (revision_v2/submission/)
- [x] revision_v2/submission/LAR_main_text_v2.docx (78348 bytes)
- [x] revision_v2/submission/LAR_additional_file_1.docx (48442 bytes)
- [x] revision_v2/submission/supplementary_tables_v2.xlsx (158803 bytes)
- [x] revision_v2/submission/response_to_reviewers_v2.docx (47057 bytes)
- [x] revision_v2/submission/cover_letter_v2.docx (38750 bytes)
- [x] figures/: 9 PDF + 9 TIFF (Fig 1–6, SFig S1–S3)
- [x] manuscript/revised_v2/abstract_structured.md (344 words)
- [x] manuscript/revised_v2/introduction.md (842 words)
- [x] manuscript/revised_v2/results.md (4839 words)
- [x] manuscript/revised_v2/discussion.md (3472 words)
- [x] manuscript/revised_v2/methods.md (4613 words)
- [x] manuscript/revised_v2/figure_legends.md (2133 words)
- [x] manuscript/revised_v2/table1.md (412 words)
- [x] manuscript/revised_v2/supplementary_table_captions.md (1457 words)
- [x] manuscript/revised_v2/supplementary_note1.md (821 words)
- [x] manuscript/revised_v2/supplementary_note2.md (1651 words)
- [x] manuscript/revised_v2/references.md (641 words)
- [x] manuscript/revised_v2/declarations.md (266 words)
- [x] manuscript/revised_v2/abbreviations.md (138 words)
- [x] manuscript/revised_v2/keywords.md (15 words)
- [x] revision_v2/text_changes.md (7767 words)
- [x] revision_v2/response_to_reviewers.md (4410 words)
- [x] revision_v2/cover_letter.md (703 words)
- [x] revision_v2/CHANGELOG.md (12459 words)

## 6. Open items for the authors (human decision)
- Appendix C, decided 2026-10-20: '+1.14' (Rps2 liver bias) and '0.005' (liver WASP − naive mean) are coincidental matches of NEW values with v1 strings (ACCEPTED contexts in appendix_c_grep.sh); 'heart' is the GEO tissue attribute of the rejected GSE55792 samples.
- The FVB_NJ_v3 comparison (optional Phase G) was not performed; the Discussion says so and names the Rps2 skin locus as the first place to look.
- github_repo working tree: the v1-era working-tree copies (seven with uncommitted author changes) were moved to revision_v2/backup_v1_worktree/ with their patches in Phase S4; committed modes restored; the v1 scripts are in legacy/v1/ on the branch.

**Result: 0 failing check(s).**
