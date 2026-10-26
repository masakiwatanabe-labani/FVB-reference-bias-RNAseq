#!/usr/bin/env python3
"""Phase S3 (2026-10-23): README.md and CITATION.cff of the public repository, generated (never typed by hand).
Numbers: «key|fmt» placeholders resolved from manuscript/revised_v2/generated/numbers.tsv (M-1 audits the output).
Dataset table: the body of manuscript/revised_v2/generated/table1.md (make_table1.py, from the sample manifest), i.e. the same
source as Table 1. Title: revision_v2/submission/title_page.md. Output: repo_release/README.md, repo_release/CITATION.cff
(copied into github_repo by commit_phase_h.sh). Usage: make_repo_readme.py
Phase S5 (2026-10-24): CITATION.cff — preferred-citation status: submitted; author e-mail removed; date-released 2026-10-24."""
import csv, re, os
P = "/usr/local/jupyter/FVB_bias"; OUT = f"{P}/repo_release"; os.makedirs(OUT, exist_ok=True)
N = {r["key"]: r["value"] for r in csv.DictReader(open(f"{P}/manuscript/revised_v2/generated/numbers.tsv"), delimiter="\t")}
def fmt(k, f):
    x = float(N[k])
    return f"{int(round(x)):,}" if f == ",d" else f"{int(round(x))}" if f == "d" else format(x, f)
def R(s): return re.sub(r"«([^|»]+)\|([^»]+)»", lambda m: fmt(m.group(1), m.group(2)), s)
title = [l.strip() for l in open(f"{P}/revision_v2/submission/title_page.md", encoding="utf-8").read().splitlines()
         if l.strip() and not l.startswith("**") and not l.startswith("#")][0]
t1 = open(f"{P}/manuscript/revised_v2/generated/table1.md", encoding="utf-8").read()
t1 = re.sub(r"^# Table 1\..*\n", "", t1, flags=re.M)
t1 = re.sub(r"^Generated from the sample manifest.*\n", "", t1, flags=re.M)
t1 = re.sub(r"^## ", "### ", t1, flags=re.M).strip()
t1 = t1.replace("(see *Sample identity verification*)", "(manuscript Methods, *Sample identity verification*)")

readme = R(f"""# FVB reference-genome mapping bias in RNA-seq

Code, sample manifest, results tables and revision records for:

> **{title}**
> Masaki Watanabe, Nobuya Sasaki. *Laboratory Animal Research* (revised manuscript).

## Versions

**The branch `v2-corrected` and the tag `v2-corrected-r1` contain the analysis reported in the revised manuscript.** The earlier
content of this repository (branch `main`) corresponds to the originally submitted version, whose dataset selection was corrected
in the revision; it is kept unchanged in the history, and its scripts are under `legacy/v1/` on this branch for reference.

## Summary

RNA-seq reads from FVB/NJ mice were quantified in three alignment arms — naive alignment to the C57BL/6J reference GRCm39, the
same alignment after WASP filtering, and alignment to an FVB/NJ-personalised GRCm39 built with g2gtools from the Mouse Genomes
Project SNPs — and the FVB/NJ-versus-C57BL/6J log₂ fold change of every gene was compared between arms.

- **Strain comparison:** «animals_fvb_samples|d» FVB/NJ and «animals_b6_samples|d» C57BL/6J samples («animals_fvb_n|d» and «animals_b6_n|d» animals) in four tissues from
  two laboratories — mesenteric adipose tissue, hypothalamus and liver from a fructose-feeding study (GSE123875, GSE123893,
  GSE123894; the three tissues come from the same animals), and back skin from an imiquimod strain panel (GSE86315).
- **F1 allele-specific expression:** «n_f1|d» F1 (C57BL/6J × FVB) samples: «animals_f1_cortex_n|d» animals from GSE293440 and «animals_f1_keratinocyte_cultures_n|d» keratinocyte cultures
  from GSE200632 (donor animals not recorded).
- **Negative control:** «n_pyvt|d» mammary tumours of the C57BL/6-congenic B6.FVB-Tg(MMTV-PyVT) line (GSE175625).
- **Catalog:** «catalog_n|d» genes whose FVB/B6 estimate depends on the reference; at «decomp_both|d» of them FVB/NJ reads are demonstrably recovered
  on the FVB/NJ reference (component `bias_matched` of the decomposition), and at the others the change mainly reflects
  C57BL/6J reads lost on that reference.
- **Survey of practice:** «survey_total|d» public GEO series with FVB RNA-seq data (of «survey_candidates_n|d» candidates, re-screened for RNA-Seq library
  strategy and FVB background per sample): «survey_stated_B6_n|d» stated a B6 reference, «survey_not_reported_n|d» did not report one,
  and none reported a strain-aware method (`survey/`).

Every animal was checked against its recorded strain by the allele fractions at FVB/NJ-homozygous exonic SNP sites before
use (sample-identity check; results in `analysis_v2/results/qc/sample_identity_verification.tsv`).

## Datasets

Generated from the sample manifest (`revision_v2/manifest/manifest.tsv`) by `revision_v2/scripts/make_table1.py`, the same
source as Table 1 of the manuscript.

{t1}

Raw reads are not included; all datasets are public in NCBI GEO / SRA (downloaded from ENA by the pipeline).

## Directory map

| Path | Contents |
|---|---|
| `revision_v2/manifest/` | sample manifest (one row per sample: accession, strain, tissue, condition, inclusion and reason) |
| `revision_v2/scripts/` | manifest construction, sample-identity check (`strain_identity_qc.py`, `A2_*`), Table 1, supplementary tables, figures (`fig*.py`, `make_figures.sh`), manuscript assembly and `final_check.sh` |
| `analysis_v2/pipeline/` | per-animal processing (`run_sample.sh`), batch control, count merging, edgeR and bias analysis, decomposition, F1 allele-specific expression, follow-up analyses |
| `analysis_v2/params.yaml` | fixed analysis parameters and reference paths |
| `analysis_v2/results/` | result tables: `qc/` (sample-identity check), `catalog/`, `predictors/`, `negative_control/`, `sensitivity_matched/`, `within_fvb/`, `f1_ase/`, `identity/`, `worked_examples/`, `survey/` |
| `analysis_v2/figures/` | data plotted in every figure panel (`*_data.tsv`) and checksums |
| `revision_v2/submission/figures/` | Figures 1–6 and Supplementary Figures S1–S3 as submitted (PDF and 600-dpi TIFF; md5 in `analysis_v2/figures/FIGURES.md5`) |
| `survey/` | GEO survey script and classification table (see `survey/README.md`) |
| `identity_check/` | the fixed 10% site subset of the sample-identity check, its traversal intervals and the script that regenerates them (see `identity_check/README.md`) |
| `scripts/pipeline/01_build_references.sh` | reference construction (GRCm39, MGP FVB/NJ variants, g2gtools personalised genome, STAR indices) |
| `manuscript/revised_v2/` | generated numbers (`generated/numbers.tsv`) and supplementary tables. The manuscript and its supplementary files are available from the journal. |
| `revision_v2/` | revision records: `CHANGELOG.md`, `text_changes.md`, `change_ledger.tsv`, audit notes (`audit/`), number freeze and independent recomputation (`freeze/`) |
| `legacy/v1/` | scripts of the originally submitted version (not used by the revised analysis) |

## How to run (pipeline order)

The scripts expect the project root in the variable `PROJ` near the top of each script; set it to your working directory.
Software versions are listed in `environment.yml` and in the manuscript Methods.

1. **References:** `scripts/pipeline/01_build_references.sh` — GRCm39 FASTA and Ensembl 110 GTF, MGP REL-2021 FVB_NJ variants,
   g2gtools VCI/patch/GTF conversion, STAR indices for both references.
2. **Manifest:** `revision_v2/scripts/A1_fetch_metadata.py`, `A1_merge_metadata.py`, `B1_fetch_fructose_samples.py`,
   `build_manifest.py` — sample attributes from GEO, ENA and BioSample, and the inclusion rules.
3. **Identity-check sites:** `revision_v2/scripts/A2_build_qc_sites.sh` (FVB/NJ-homozygous exonic SNP sites; verdict rules in
   `strain_identity_qc.py`), then `identity_check/make_identity_subset.py` for the fixed 10% subset used for every animal except
   the F1 hybrids.
4. **Per-animal processing:** `analysis_v2/pipeline/run_batch.sh <batch>` (or `scheduler.sh`), which runs `run_sample.sh` per animal:
   ENA download and md5 → fastp → STAR naive arm with WASP tags → sample-identity check (stops on a mismatch) → WASP filter →
   STAR FVB/NJ arm → featureCounts (unique-only and multi-mapping).
5. **Gene metadata:** `analysis_v2/pipeline/gene_meta.sh`.
6. **Primary analysis:** `analysis_v2/pipeline/final_run.sh` — `merge_counts.py` → `bias_analysis.py` (edgeR via `run_edger.R`;
   catalog, predictors, within-FVB comparisons) → `revision_v2/scripts/make_figures.sh`.
7. **Further analyses** in `analysis_v2/pipeline/`: `pyvt_fullsite_recount.sh` → `followup_analysis.py` (negative control,
   matched-reference sensitivity) → `bias_decomposition.py`; `catalog_clusters_exonic_snps.py`, `rps2_*.py`.
8. **F1 allele-specific expression:** `f1_ase_rerun.sh` → `f1_ase_summary_v2.py`, `f1_wasp_allele_removal.py`,
   `f1_sites_catalog_vs_other.py`.
9. **Survey:** `survey/rescreen_fetch.py` (GEO/SRA records cached in `survey/raw/`), `survey/rescreen_classify.py` (re-screen and
   classification with `classify()` of `survey/stage3_survey_v2.py`); Fig. 1 by `revision_v2/scripts/fig1_survey.py`.
10. **Numbers, tables and checks:** `revision_v2/scripts/results_numbers.py` (`numbers.tsv`), `make_table1.py`,
    `build_supplementary_tables.py`, then `build_submission.sh build` and `final_check.sh`.

## Citation

Please cite the article (see `CITATION.cff`):
Watanabe M, Sasaki N. {title}. *Laboratory Animal Research* (in revision).

## License

MIT (see `LICENSE`).
""")
open(f"{OUT}/README.md", "w", encoding="utf-8").write(readme)

cff = f"""cff-version: 1.2.0
message: "If you use this software or the bias catalog, please cite the article below."
type: software
title: "FVB reference-genome mapping bias in RNA-seq: analysis code"
version: "v2-corrected-r1"
abstract: >
  Code, sample manifest and result tables for the revised manuscript: three-arm RNA-seq quantification (naive GRCm39,
  WASP-filtered, FVB/NJ-personalised reference) of FVB/NJ and C57BL/6J animals, F1 allele-specific expression, a
  negative control, the gene catalog and its decomposition, and a survey of reference-genome practice.
authors:
  - family-names: "Watanabe"
    given-names: "Masaki"
  - family-names: "Sasaki"
    given-names: "Nobuya"
repository-code: "https://github.com/masakiwatanabe-labani/FVB-reference-bias-RNAseq"
license: MIT
date-released: "2026-10-24"
keywords:
  - RNA-seq
  - reference bias
  - FVB mouse
  - personalised reference genome
  - g2gtools
  - WASP
  - featureCounts
  - SNP-mediated bias
preferred-citation:
  type: article
  title: "{title}"
  authors:
    - family-names: "Watanabe"
      given-names: "Masaki"
    - family-names: "Sasaki"
      given-names: "Nobuya"
  journal: "Laboratory Animal Research"
  year: 2026
  status: submitted
"""
open(f"{OUT}/CITATION.cff", "w", encoding="utf-8").write(cff)
print("README.md", len(readme.split()), "words; CITATION.cff written")
