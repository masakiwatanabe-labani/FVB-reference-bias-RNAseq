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
open(f"{F}/FREEZE.md", "w").write("\n".join(L)); print("\n".join(L[:12]))
