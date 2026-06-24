# FVB Reference-Genome Mapping Bias in RNA-seq

Analysis code and figure-generation scripts for:

> **"SNP-mediated reference-genome mapping bias in FVB RNA-seq: a systematic characterisation
> across five tissue datasets with a 32-gene bias catalog"**
> *(manuscript in preparation)*

---

## Overview

This repository contains the analysis pipeline and figure scripts used to characterise
SNP-mediated reference-genome mapping bias in FVB versus C57BL/6J (B6) RNA-seq experiments.
Raw FASTQ/BAM files are **not** included; all source datasets are publicly available on
NCBI GEO (accession numbers below). The code is provided for transparency and reproducibility.

---

## GEO datasets

| GEO accession | Tissue | Substrain | Library | n FVB | n B6 |
|---|---|---|---|---|---|
| GSE123875 | Mesenteric adipose | FVB/NJ | Paired-end | 6 | 6 |
| GSE123893 | Liver | FVB/NJ | Paired-end | 6 | 6 |
| GSE123894 | Kidney | FVB/NJ | Paired-end | 6 | 6 |
| GSE135230 | Cardiac muscle | FVB/N | Single-end | 10 | 12 |
| GSE175625 | Mammary tumour (PyVT) | FVB/N | Single-end | 7 | 6 |
| GSE293440 | Frontal cortex (F1) | FVB/NJ | Paired-end | 8 | — |
| GSE200632 | Keratinocytes (F1) | FVB/N | Paired-end | 4 | — |

Reference genome: GRCm39 (Ensembl release 110).  
FVB/NJ variant calls: Mouse Genomes Project REL-2021 (`REL-2112-v8-SNPs_Indels/`).

---

## Software versions

| Tool | Version | Purpose |
|---|---|---|
| STAR | 2.7.10a | Read alignment (B6 and FVB reference) |
| fastp | 1.3.3 | Read quality control and adapter trimming |
| samtools | 1.17 | BAM sorting, filtering, indexing |
| featureCounts (Subread) | 2.0.3 | Gene-level read counting (unique-only primary; −M −−fraction for sensitivity) |
| g2gtools | 0.2.7 | FVB/NJ personalised reference construction |
| GATK ASEReadCounter | 4.6.2.0 | Allele-specific expression at SNP sites |
| Python | 3.11 | Analysis and figure scripts |
| NumPy | 2.4.6 | |
| SciPy | 1.17.1 | |
| pandas | 3.0.3 | |
| matplotlib | 3.11.0 | |
| pysam | 0.22 | VCF/BAM access from Python |
| R | 4.2.2 | edgeR differential expression |
| edgeR | 3.40.2 | TMM normalisation and quasi-likelihood F-tests |

---

## Repository structure

```
github_repo/
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml          # conda environment (Python dependencies)
├── scripts/
│   ├── pipeline/
│   │   ├── 01_build_references.sh       # VCF filtering → g2gtools VCI → STAR indexes
│   │   ├── 02_align_wasp_featurecounts.sh  # 3-arm alignment + featureCounts
│   │   ├── 03_f1_ase.sh                 # F1 ASE alignment + GATK ASEReadCounter
│   │   ├── 04_bias_analysis.py          # CPM normalisation → bias table → catalog
│   │   └── 05_multimapping_sensitivity.py  # −M −−fraction re-counting + Δ/mechanism
│   └── figures/
│       ├── fig01_survey.py
│       ├── fig02_f1_ase.py
│       ├── fig03_strain_comp.py
│       ├── fig04_predictors.py
│       ├── fig05_catalog.py
│       ├── fig06_worked_example.py
│       └── sfig_s1_multimap.py
└── tables/
    └── README_tables.md     # Column descriptions for Supplementary Tables S1–S4
```

---

## Analysis execution order

Set the project root before running any script:

```bash
export FVB_PROJ=/path/to/your/working/directory
```

Then run the pipeline steps in order:

```bash
# 1. Build reference genome and FVB/NJ personalised reference
bash scripts/pipeline/01_build_references.sh

# 2. Three-arm alignment (naive / WASP / FVB reference) + featureCounts
#    Download FASTQs from SRA before this step (see SRR accessions in script header)
bash scripts/pipeline/02_align_wasp_featurecounts.sh

# 3. F1 ASE analysis (optional; requires GSE293440 and GSE200632 FASTQs)
bash scripts/pipeline/03_f1_ase.sh

# 4. Bias quantification and catalog construction (outputs to $FVB_PROJ/results/stage6/)
python scripts/pipeline/04_bias_analysis.py

# 5. Multi-mapping sensitivity analysis (Supplementary Note 2)
python scripts/pipeline/05_multimapping_sensitivity.py
```

### Regenerating figures

```bash
# Main figures (requires step 4 to have completed)
python scripts/figures/fig01_survey.py
python scripts/figures/fig02_f1_ase.py
python scripts/figures/fig03_strain_comp.py
python scripts/figures/fig04_predictors.py
python scripts/figures/fig05_catalog.py
python scripts/figures/fig06_worked_example.py

# Supplementary Figure S1 (requires step 5 to have completed)
python scripts/figures/sfig_s1_multimap.py
```

All figures are saved as PNG (200 dpi) and PDF to `$FVB_PROJ/results/figures/`.

---

## Important notes

**Raw data**: FASTQ files must be downloaded from SRA/ENA before running the alignment step.
The SRR accession lists for each dataset are at the top of
`scripts/pipeline/02_align_wasp_featurecounts.sh`. Tools such as `prefetch` + `fasterq-dump`
(SRA Toolkit) or `aria2c` with ENA FTP URLs can be used.

**Disk usage**: Peak disk use during alignment is approximately 150 GB (raw FASTQs + BAMs
for one dataset). BAMs and trimmed FASTQs should be deleted after featureCounts runs per
dataset to manage disk space.

**g2gtools Python 3 compatibility**: g2gtools v0.2.7 requires a two-line patch to
`chain.py` to replace deprecated `.next()` iterator calls with Python 3 `next()`.
See `scripts/pipeline/01_build_references.sh` for details.

**PyVT tumour dataset**: GSE175625 (mammary, FVB/N) carries somatic copy number variants
driven by the MMTV-PyVT transgene. Bias values for this dataset are retained as a separate
annotation column (`pvt_bias`) in Supplementary Table S2 but are excluded from the primary
catalog criterion (see Methods).

**FVB/N vs FVB/NJ substrain mismatch**: The personalised reference was built from FVB/NJ
(JAX) variant calls. GSE135230 and GSE175625 used FVB/N mice; FVB/N-private variants are
not corrected by this VCI. See Methods and Discussion for details.

**Manuscript version**: The version of this repository corresponding to the submitted
manuscript is archived on Zenodo at [DOI — to be added upon acceptance].

---

## Citation

If you use these scripts or the bias catalog, please cite the manuscript (citation below)
and this repository. A machine-readable citation is provided in `CITATION.cff`.

---

## License

Code: MIT License (see `LICENSE`).  
Bias catalog (Supplementary Table S2): CC BY 4.0.
