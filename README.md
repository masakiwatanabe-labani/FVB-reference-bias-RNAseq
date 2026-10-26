# FVB reference-genome mapping bias in RNA-seq

Code, sample manifest, results tables and revision records for:

> **Reproducible gene-level reference-genome mapping bias in FVB mouse RNA-seq: implications for animal study reproducibility**
> Masaki Watanabe, Nobuya Sasaki. *Laboratory Animal Research* (revised manuscript).

## Versions

**The branch `v2-corrected` and the tag `v2-corrected-r1` contain the analysis reported in the revised manuscript.** The earlier
content of this repository (branch `main`) corresponds to the originally submitted version, whose dataset selection was corrected
in the revision; it is kept unchanged in the history, and its scripts are under `legacy/v1/` on this branch for reference.

## Summary

RNA-seq reads from FVB/NJ mice were quantified in three alignment arms — naive alignment to the C57BL/6J reference GRCm39, the
same alignment after WASP filtering, and alignment to an FVB/NJ-personalised GRCm39 built with g2gtools from the Mouse Genomes
Project SNPs — and the FVB/NJ-versus-C57BL/6J log₂ fold change of every gene was compared between arms.

- **Strain comparison:** 36 FVB/NJ and 34 C57BL/6J samples (23 and 23 animals) in four tissues from
  two laboratories — mesenteric adipose tissue, hypothalamus and liver from a fructose-feeding study (GSE123875, GSE123893,
  GSE123894; the three tissues come from the same animals), and back skin from an imiquimod strain panel (GSE86315).
- **F1 allele-specific expression:** 12 F1 (C57BL/6J × FVB) samples: 8 animals from GSE293440 and 4 keratinocyte cultures
  from GSE200632 (donor animals not recorded).
- **Negative control:** 13 mammary tumours of the C57BL/6-congenic B6.FVB-Tg(MMTV-PyVT) line (GSE175625).
- **Catalog:** 138 genes whose FVB/B6 estimate depends on the reference; at 57 of them FVB/NJ reads are demonstrably recovered
  on the FVB/NJ reference (component `bias_matched` of the decomposition), and at the others the change mainly reflects
  C57BL/6J reads lost on that reference.
- **Survey of practice:** 656 public GEO series with FVB RNA-seq data (of 1369 candidates from a broad GEO search, screened
  for RNA-Seq library strategy and FVB background per sample): 620 stated a B6 reference, 30 did not report one,
  and 6 reported a strain-aware method (`survey/`; the series deposited by the authors' group is not counted).

Every animal was checked against its recorded strain by the allele fractions at FVB/NJ-homozygous exonic SNP sites before
use (sample-identity check; results in `analysis_v2/results/qc/sample_identity_verification.tsv`).

## Datasets

Generated from the sample manifest (`revision_v2/manifest/manifest.tsv`) by `revision_v2/scripts/make_table1.py`, the same
source as Table 1 of the manuscript.

### Three-arm strain comparison datasets (FVB/NJ versus C57BL/6J)

| GEO accession | Tissue | Cohort | *n* FVB/NJ (conditions) | *n* C57BL/6J (conditions) | Library | Read length | Runs per sample | Sequence identity QC passed |
|---|---|---|---|---|---|---|---|---|
| GSE123875 | Mesenteric adipose | Fructose study (PRJNA510125; SubSeries of GSE123896) | 8 (fructose ×5; water ×3) | 6 (fructose ×2; water ×4) | Paired-end | 69 bp | 3 | 14/14 |
| GSE123893 | Hypothalamus | Fructose study (PRJNA510129; SubSeries of GSE123896) | 8 (fructose ×4; water ×4) | 8 (fructose ×4; water ×4) | Paired-end | 101 bp | 3 | 16/16 |
| GSE123894 | Liver | Fructose study (PRJNA510133; SubSeries of GSE123896) | 12 (fructose ×6; water ×6) | 12 (fructose ×6; water ×6) | Paired-end | 69 bp | 3 | 24/24 |
| GSE86315 | Back skin | Imiquimod strain panel (PRJNA341396) | 8 (IMQ F ×2; IMQ M ×2; control F ×2; control M ×2) | 8 (IMQ F ×2; IMQ M ×2; control F ×2; control M ×2) | Single-end | 52 bp | 1 | 16/16 |

### Negative-control dataset (no FVB/NJ-versus-C57BL/6J contrast)

| GEO accession | Tissue | Strain (GEO attribute) | *n* | Conditions | Library | Read length | Sequence identity QC |
|---|---|---|---|---|---|---|---|
| GSE175625 | Mammary tumour (MMTV-PyVT) | B6.FVB-Tg(MMTV-PyVT)634Mul/LellJ (C57BL/6 congenic) | 13 | Methylparaben ×4; Propylparaben ×6; Vehicle ×3 | Single-end | 101 bp | 13/13 (expected: C57BL/6) |

### F1 hybrid datasets used for allele-specific expression analysis (Fig. 2)

| GEO accession | Tissue | Cross (GEO attribute) | *n* (F1 samples) | Library | Read length | Sequence identity QC |
|---|---|---|---|---|---|---|
| GSE293440 | Cerebral (frontal) cortex | B6xFVB | 8 | Paired-end | 150 bp | 8/8 (expected: F1) |
| GSE200632 | Primary keratinocytes (LSL-KrasG12D, Cre-transduced) | C57BL/6:FVB/N (F1:50/50) | 4 | Paired-end | 151 bp | 4/4 (expected: F1) |

Adipose, hypothalamus and liver samples of the fructose study were collected from the same animals; *n* per tissue counts samples (fructose study: FVB/NJ 28 samples from 15 animals, C57BL/6J 26 samples from 15 animals; all strain-comparison samples: FVB/NJ 36 samples from 23 animals, C57BL/6J 34 samples from 23 animals; Supplementary Table S1). The PyVT tumours are one primary tumour per mouse; GSE200632 *n* counts keratinocyte cultures from newborn pups (number of donor animals not recorded).
Sequence identity QC: allele fractions at FVB/NJ-homozygous exonic SNP sites (203,510; a fixed 10% subset for all but the F1 samples) in the B6-reference alignment (manuscript Methods, *Sample identity verification*); full per-sample results in Supplementary Table S1.
Sex and age are not recorded in the GEO/BioSample attributes of the fructose study; the imiquimod study records sex (M/F). All FVB/NJ animals are from the Jackson Laboratory substrain, matched to the FVB_NJ variant set.

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
9. **Survey:** `survey/broad_search.py` (candidate list), `survey/rescreen_fetch.py` (GEO/SRA records cached in `survey/raw/`), `survey/rescreen_classify.py` (re-screen and
   classification with `classify()` of `survey/stage3_survey_v2.py`); Fig. 1 by `revision_v2/scripts/fig1_survey.py`.
10. **Numbers, tables and checks:** `revision_v2/scripts/results_numbers.py` (`numbers.tsv`), `make_table1.py`,
    `build_supplementary_tables.py`, then `build_submission.sh build` and `final_check.sh`.

## Citation

Please cite the article (see `CITATION.cff`):
Watanabe M, Sasaki N. Reproducible gene-level reference-genome mapping bias in FVB mouse RNA-seq: implications for animal study reproducibility. *Laboratory Animal Research* (in revision).

## License

MIT (see `LICENSE`).
