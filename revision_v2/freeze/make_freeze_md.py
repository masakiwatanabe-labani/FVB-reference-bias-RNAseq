#!/usr/bin/env python3
"""Phase M-4: write revision_v2/freeze/FREEZE.md from the freeze records (token_audit.tsv, recompute_report.tsv, inputs.sha256,
deliverables.sha256) and the before-hashes of the authorised exceptions (before_hashes.tsv, recorded from the pre-change files).
Usage: make_freeze_md.py <parent commit>"""
import hashlib, os, sys, datetime
import pandas as pd
P = "/usr/local/jupyter/FVB_bias"; F = f"{P}/revision_v2/freeze"
def sha(p): return hashlib.sha256(open(f"{P}/{p}", "rb").read()).hexdigest() if os.path.exists(f"{P}/{p}") else "(absent)"
tok = pd.read_csv(f"{F}/token_audit.tsv", sep="\t", dtype=str, keep_default_na=False)
rep = pd.read_csv(f"{F}/recompute_report.tsv", sep="\t", dtype=str, keep_default_na=False)
bef = pd.read_csv(f"{F}/before_hashes.tsv", sep="\t", dtype=str, keep_default_na=False)
parent = sys.argv[1] if len(sys.argv) > 1 else "?"
L = [f"# Number freeze — Phase M ({datetime.date.today().isoformat()})", "",
     f"Sealed in the commit that adds this file (parent {parent}, branch v2-corrected); that commit is the one the authors tag v2-corrected-r1.", "",
     "## What was verified", "",
     "- **M-0 inputs frozen**: `inputs.sha256` lists every frozen input (numbers.tsv, supplementary_tables_v2.xlsx, every file of RESULTS.md5 and FIGURES.md5, and the two md5 lists); final_check fails if any changes. `build_submission.sh` does not rewrite the xlsx or the figures once frozen.",
     f"- **M-1 numeric-token audit** (`token_audit.py`, `resolve_tokens.py` → `token_audit.tsv`): {len(tok)} numeric tokens in the main text (all sections, Table 1, legends), Additional file 1, both letters and the text of the nine figure PDFs.",
     "- **M-2 independent recomputation** (`recompute_independent.py` → `recompute_report.tsv`): no import of the pipeline or of the text/number generators, no numbers.tsv; primary result files only (edgeR tables, count matrices and featureCounts summaries, STAR junction lists and Log.out, F1 and identity per-site tables, WASP allele-removal table, within-strain DE tables, FVB_NJ VCFs and the Ensembl 110 GTF via bcftools/bedtools, rps2_* files, vW tag counts, survey table, submitted-version catalog for the letters). Bounds pass only if the exact value satisfies the bound as printed.",
     "- **D3 double-rounding scan** (`double_rounding_scan.py`): every key printed in the text compared with the value rounded once from full precision (`generated/numbers_fullprecision.tsv`).",
     "- **Seal**: `deliverables.sha256` (every file in revision_v2/submission/, figures included).", "",
     "## Token classes (M-1)", "", "| class | tokens |", "|---|---|"]
L += [f"| {k} | {v} |" for k, v in tok["class"].value_counts().items()]
L += ["", "KEY = equals a numbers.tsv value at the printed precision; XLSX = derivable from the supplementary tables; STRUCT = citation, label, accession, year, version, parameter or fixed threshold; M2 = not in numbers.tsv or the tables but reproduced by a PASS claim of the independent recomputation (claim id in `ref`); UNMAPPED = see below.", "",
      "UNMAPPED (2): **1,047** (GEO search hits, Methods) — carried from the submitted version; no search log exists in the project; not recomputable. Plausibility check 2026-10-23 through E-utilities (db=gds, `FVB AND \"expression profiling by high throughput sequencing\"[DataSet Type] AND \"Mus musculus\"[Organism] AND gse[Entry Type]`): 1,073 series today (≥ 1,047, as expected for a later date); text unchanged. **290** (size of the GSE293440 series, Methods) — verified against the GEO SOFT record (290 `!Series_sample_id` entries, 2026-10-23).", "",
      "## Claims (M-2)", "", "| group | PASS | FAIL | not recomputable |", "|---|---|---|---|"]
for g_, d in rep.groupby("group"):
    L.append(f"| {g_} | {(d.status == 'PASS').sum()} | {(d.status == 'FAIL').sum()} | {(d.status == 'NOT_RECOMPUTABLE').sum()} |")
L.append(f"| **total** | **{(rep.status == 'PASS').sum()}** | **{(rep.status == 'FAIL').sum()}** | **{(rep.status == 'NOT_RECOMPUTABLE').sum()}** |")
L += ["", "## Authorised exceptions to the freeze (authors' decisions, 2026-10-23)", "",
      "| id | change | file | sha256 before | sha256 after |", "|---|---|---|---|---|"]
for r in bef.itertuples():
    L.append(f"| {r.id} | {r.change} | `{r.file}` | `{r.before[:16]}…` | `{sha(r.file)[:16]}…` |")
L += ["", "Full before/after hashes: `before_hashes.tsv` (before) and `inputs.sha256` / `deliverables.sha256` (after). Every other file of FIGURES.md5 and RESULTS.md5 is byte-identical to the freeze of 2335065; `rps2_pseudogene_alleles.tsv` is unchanged (no SNP inside the copied part is new).", "",
      "- **D1** Fig. 1 regenerated (PDF, TIFF; `_data.tsv` byte-identical): `fig1_survey.py` read the submitted-version F1 table; panel B now prints 0.480 / 0.512 (pdftotext), no 0.333 / 0.386.",
      "- **D1b** Supplementary Fig. S1 regenerated: reference animals by rule (most sites of depth ≥ 20 per strain among the v2 animals evaluated at the full site set; candidates in `analysis_v2/results/identity/full_site_tables/candidates.tsv`): adipose_B6_CWM8 (35,058 sites) and adipose_FVB_FWM9 (30,685), both GSE123875; the six per-site tables (3.1 MB in all) are in that directory, the only input of `sfig_gse55792_blocks.py` besides the reference .fai.",
      "- **D2** Rps2 exonic variants recounted with exonic = union of the exons of all Ensembl 110 Rps2 transcripts (`analysis_v2/results/worked_examples/rps2_exonic_variants_union.tsv`): 16 SNPs, 5 indels; 5 SNPs inside the part copied by the pseudogenes (ps10 5/5, ps13 5/5, ps6 3/5 FVB/NJ), 11 outside; 0 indels inside (case 1). The genomic blastn cross-check differs only at 17:24,940,643 (a 106-bp exon below e < 1e-50 as a standalone HSP; inside the copied part in the mature-sequence alignment).",
      "- **D3** Rps2/pseudogene CPM keys stored at full precision; Note 2 table 209.75 → 209.76. Scan: the only double-rounding defect before the change was this value; after: 0.",
      "- **D4 / E2 / M-3a–c / D5 / E1** text changes (trigger PHASE-M; text_changes.md, Phase M). New numbers.tsv keys: `vcf_fvbnj_homalt_autosomal_indels`, `identity_subset_validation_n_animals` (14) and `_max_absdiff` (0.0135), `nc_catalog_genes_expressed_min/max` (76/86), `decomp_sum_mismatch_n_genes` (1) and `decomp_Gm26740_sum_minus_mean_bias` (−0.0186), `sfig1_shared_block_chromosomes` (2, 3, 8, 10 and 11), `sfig1_shared_block_n_chromosomes`, `sfig1_ref_b6_max_hom_fvb` (0.0326), `sfig1_ref_fvb_min_hom_fvb` (0.691). Changed values: `rps2_snps_exonic` 13 → 16, `rps2_indels_exonic` 4 → 5, `rps2_exonic_snps_outside_pseudogene_sequence` 8 → 11; 59 CPM keys changed in stored precision only.",
      "- **D6** `manuscript/revised_v2/generated/wasp_attrition_by_group.tsv` (input of Sheet S4b): input outside analysis_v2/results/, v2-derived, independently verified (M-2 group I, 9/9 PASS); tracked in git; the xlsx was not rebuilt.", "",
      "## Other records", "",
      "- Duplicate keys in numbers.tsv (`ho_adipose/hypothalamus/liver/skin_n_catalog`, three rows each; `pred_biotype_kruskal_p`, two rows): identical values in all rows (E3); left unchanged during the freeze.",
      "- D1.1 grep after the changes: no figure or table script reads outside analysis_v2/results/, revision_v2/manifest/ and ref/ except Sheet S4b (D6).",
      "- `revision_v2/freeze/cache/` (intermediate BEDs of M-2, 431 MB) is not committed.",
      "- `submission_checklist.md` was generated before the seal; its seal line therefore reads \"not yet written\"; the final_check run after the seal (final_check.md) verifies it.", ""]
if os.path.exists(f"{F}/before_hashes_N.tsv"):
    bn = pd.read_csv(f"{F}/before_hashes_N.tsv", sep="\t", dtype=str, keep_default_na=False)
    sz = pd.read_csv(f"{F}/figure_sizes.tsv", sep="\t", dtype=str)
    L += ["## Phase N (2026-10-23): fixes from the package check — authorised exceptions", "",
          "No number changed: every _data.tsv is byte-identical to the frozen version; the M-1 token audit of the figure PDFs differs only by layout tokens (axis ticks, digits of rotated labels split by pdftotext), the Fig. 3A key (117 of 138 catalog genes tested in liver; n = 40/31/46 tested in liver) and Fig. 4 p values printed as \"p < 1e-15\" (was 0.0e+00, 2.7e-17, 1.6e-25, 1.8e-35); all are claims of the M-2 report (group N).", "",
          "- N-1a Fig. 6D title: \"largest single-tissue bias of the direct-misalignment class (skin +5.36 log₂)\"; Supplementary Fig. S2D title: \"two largest |Δ| of the catalog\" (the only other \"largest\" in a figure).",
          "- N-1b Fig. 4: printed p below 1e-15 shown as \"p < 1e-15\" (no other figure printed a p of 0).",
          "- N-1c Supplementary Fig. S1 legend: \"(both mesenteric adipose tissue, GSE123875; …)\"; bin filter stated.",
          "- N-1d Cover letter and response letter: \"(Fig. 1A) was not affected; the schematic in Fig. 1B, which quoted the F1 values, has been updated\".",
          "- N-2a Results and Note 2: \"the other eleven exonic SNPs lie outside the part that the pseudogenes copy\" (10 lie in exon segments absent from ENSMUST00000146867, the transcript the pseudogenes copy; 17:24,940,701 lies in that transcript just past the copied span).",
          "- N-2b Fig. 3A key: \"(n = … tested in liver)\" and \"catalog genes tested in liver: 117 of 138\".",
          "- N-2c TIFFs written as RGB on white (no alpha).",
          "- N-3 every figure re-laid out at its print size, all text ≥ 7 pt (sizes measured from the figure at save time):", "",
          "| figure | width (mm) | height (mm) | smallest text (pt) |", "|---|---|---|---|"]
    L += [f"| {r.figure} | {r.width_mm} | {r.height_mm} | {r.min_text_pt} |" for r in sz.itertuples()]
    L += ["", "| id | change | file | sha256 before | sha256 after |", "|---|---|---|---|---|"]
    L += [f"| {r.id} | {r.change} | `{r.file}` | `{r.before[:16]}…` | `{sha(r.file)[:16]}…` |" for r in bn.itertuples()]
    L += [""]
if os.path.exists(f"{F}/before_hashes_O.tsv"):
    bo = pd.read_csv(f"{F}/before_hashes_O.tsv", sep="\t", dtype=str, keep_default_na=False)
    oc = pd.read_csv(f"{F}/figure_overlap_check.tsv", sep="\t", dtype=str, keep_default_na=False)
    sz = pd.read_csv(f"{F}/figure_sizes.tsv", sep="\t", dtype=str)
    L += ["## Phase O (2026-10-23): figure label and occlusion fixes — authorised exceptions", "",
          "Figure-only changes; no data, number or text outside the figures changed: the letters, the xlsx and the .md/.plain.txt sources of both main docx files are byte-identical to Phase N; the two main docx files were re-packed by pandoc during the run (same text, different zip metadata). Every _data.tsv byte-identical; M-1 token differences in the figure PDFs are layout tokens only (axis ticks; Fig. 4D counts now extracted whole); M-2 368 PASS, 0 FAIL.", "",
          "- O-1a Fig. 3A: only the worked-example genes Ndufs5 and Rps2 labelled, leader lines into the empty upper left, white text background.",
          "- O-1b Fig. 5B: labels stacked beside their column with leader lines, all inside the axes.",
          "- O-1c Supplementary Fig. S2B: the ten labels stacked right of the points with leader lines.",
          "- O-1d Supplementary Fig. S2C: the residual-fraction lists moved right of the bar area (no text on any bar).",
          "- O-1e Supplementary Fig. S2D: horizontal single-letter tissue ticks (A, H, L, S, explained in the legend); one slot per tissue so that the unique and multi-mapping bars of neighbouring tissues no longer overlap.",
          "- O-2a Fig. 2C legend to the upper right (clear of the peak and the mean lines); O-2b Fig. 4A Ndufs5 label in the empty space at x < 0 with a leader line; O-2c Fig. 3B 95th-percentile triangles larger with a black edge; O-2d Fig. 1B solid \"maps\" arrow from the blue B6-allele box, above the orange box.",
          "- Also found by the automatic check and fixed: horizontal two-line tissue ticks in Figs. 3, 5A and 6 (rotated labels overlapped); Fig. 4D tick labels vertical; Fig. 2 panel letter clear of the tick labels; S2A summary box and legend clear of each other and of the bars.",
          "- O-3 `fig_common.overlap_check` (called by `save()`; a figure is not written if it fails): (i) no two visible text artists overlap (rendered text extents, leader lines excluded); (ii) no text with a background box and no legend frame overlaps a bar or patch of its axes. Exceptions: colour-bar labels and tick labels; texts inside one legend. Every figure was also rendered at 150 dpi and inspected.", "",
          "| figure | width (mm) | height (mm) | smallest text (pt) | texts checked | overlaps |", "|---|---|---|---|---|---|"]
    m = oc.set_index("figure")
    L += [f"| {r.figure} | {r.width_mm} | {r.height_mm} | {r.min_text_pt} | {m.texts_checked.get(r.figure,'')} | {m.overlaps.get(r.figure,'')} |" for r in sz.itertuples()]
    L += ["", "| id | change | file | sha256 before | sha256 after |", "|---|---|---|---|---|"]
    L += [f"| {r.id} | {r.change} | `{r.file}` | `{r.before[:16]}…` | `{sha(r.file)[:16]}…` |" for r in bo.itertuples()]
    L += [""]
if os.path.exists(f"{F}/before_hashes_P.tsv"):
    bp = pd.read_csv(f"{F}/before_hashes_P.tsv", sep="\t", dtype=str, keep_default_na=False)
    L += ["## Phase P (2026-10-23): answering the reviewers more directly — authorised exceptions", "",
          "New numbers only through numbers.tsv and M-2 (group P): share of tested genes without an exonic FVB/NJ SNP (48.1%; exons = union of all Ensembl 110 exons of the gene, SNPs = every FVB_NJ SNP record); catalog clusters by a fixed rule (single linkage of catalog genes on one chromosome with Ensembl gene starts < 1 Mb apart, clusters of ≥ 3): 39 genes in 10 clusters, 26 of them in the SNP-induced multi-mapping class; MHC 9 genes in two clusters, replication-dependent histones 11 in two clusters, Mup 3, Ifi200 3. Literature values (30,048 structural variants ≥ 50 bp; the Pde6b nonsense mutation and the 2-bp Hc deletion; ≈50× coverage on mm9) checked against the verbatim sentences of Wong et al. 2012 (revision_v2/audit/wong2012_excerpts.txt). FVB_NJ_v3 checked with NCBI Datasets: GCA_921998635.2, PacBio, falcon-kit 1.8.1, chromosome level, EMBL-EBI, released 2022-07-01, BioProject PRJEB47108 (no publication linked in NCBI; the BioProject is the data deposit of Helmy et al. 2026, Cell Genom 6:101074, whose assemblies are named {strain}_v3). Supplementary Table S2 gained the column cluster; every other cell of every sheet is identical (cell-by-cell comparison). Letters: every quotation of the revised text is checked verbatim against the main text (final_check). P-6 (FVB_NJ_v3 quantification) not run (no approval).", "",
          "| id | change | file | sha256 before | sha256 after |", "|---|---|---|---|---|"]
    L += [f"| {r.id} | {r.change} | `{r.file}` | `{r.before[:16]}…` | `{sha(r.file)[:16]}…` |" for r in bp.itertuples()]
    L += [""]
if os.path.exists(f"{F}/before_hashes_Q.tsv"):
    bq = pd.read_csv(f"{F}/before_hashes_Q.tsv", sep="\t", dtype=str, keep_default_na=False)
    L += ["## Phase Q (2026-10-23): claims aligned with the decomposition; FVB_NJ_v3 rationale; within-strain wording; Table 1 — authorised exceptions", "",
          "Text, letters and the Table 1 layout only; numbers.tsv, the xlsx and every figure unchanged. New text numbers (57 of 138; Ndufs5 bias_matched 1.47–4.61 and b6_loss +1.85 to +5.63; 598 and 539 DE genes; the v1 F1 values in the cover letter) are M-2 claims of group Q. Ensembl annotation of FVB_NJ_v3 verified online (Ensembl release 115, Mus_musculus_fvbnj.FVB_NJ_v3.115.gtf.gz; genebuild last updated 2023-06; own gene identifiers ENSMUSG00205…); cited as reference 'Ensembl FVB 2025'. Table 1: pandoc's cell paragraph style 'Compact' is absent from the localised reference docx, which made LibreOffice collapse the table; cell paragraphs now use the default style, every table has explicit dxa widths, the Table Grid style, 8-pt text and a repeated header row, and Table 1 is in its own landscape section; final_check converts the docx with LibreOffice and checks the layout (freeze/table1_libreoffice_check.png).", "",
          "| id | change | file | sha256 before | sha256 after |", "|---|---|---|---|---|"]
    L += [f"| {r.id} | {r.change} | \`{r.file}\` | \`{r.before[:16]}…\` | \`{sha(r.file)[:16]}…\` |" for r in bq.itertuples()]
    L += [""]
open(f"{F}/FREEZE.md", "w").write("\n".join(L)); print("\n".join(L[:12]))
