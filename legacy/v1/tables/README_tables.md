# Supplementary Tables

Supplementary tables are provided as Excel (.xlsx) files alongside the manuscript.
Column descriptions are listed below for each table.

---

## Supplementary Table S1 — GEO survey of reference genome usage

One row per GEO series (n = 154).

| Column | Description |
|--------|-------------|
| `gse` | GEO series accession |
| `title` | Series title (from SOFT metadata) |
| `category` | Classification: `B6_stated`, `not_stated`, or `strain_aware` |
| `data_processing_excerpt` | Relevant text from `!Sample_data_processing` or `!Series_data_processing` |

---

## Supplementary Table S2 — Bias catalog (32 genes)

One row per catalog gene. The primary catalog criterion is bias ≥ 0.5 log₂ units in ≥ 2 of
the four non-tumour tissue datasets (adipose, liver, kidney, cardiac).

| Column | Description |
|--------|-------------|
| `gene_id` | Ensembl gene ID (GRCm39 / Ensembl 110) |
| `gene_name` | Gene symbol |
| `biotype` | Ensembl gene biotype |
| `n_datasets_4t` | Number of non-tumour datasets meeting bias ≥ 0.5 threshold (0–4) |
| `n_fvbnj_support` | Number of FVB/NJ-substrain-matched datasets meeting threshold (adipose, liver, kidney) |
| `mean_bias` | Mean bias across qualifying non-tumour tissues (log₂FC_fvb − log₂FC_naive) |
| `pvt_bias` | Bias in GSE175625 (PyVT mammary tumour; may reflect CNV effects) |
| `pvt_discordant` | `True` if four-tissue mean bias is positive but PyVT bias is negative |
| `tissues_4t` | Comma-separated list of non-tumour tissues meeting threshold |
| `bias_multimap` | Mean bias under multi-mapping-inclusive counting (featureCounts -M --fraction) |
| `delta_multimap_minus_unique` | Δ = bias_multimap − mean_bias; negative = SNP-induced multi-mapping mechanism |
| `residual_fraction` | bias_multimap / mean_bias; near 1 = direct misalignment; near 0 = multi-mapping |
| `mechanism_class` | `SNP-induced multi-mapping` (|Δ| > 0.3 log₂) or `Direct SNP misalignment` |
| `mechanism_note` | Brief mechanistic interpretation |

---

## Supplementary Table S3 — Catalog sensitivity to threshold choice

One row per threshold combination. Columns record the number of catalog genes identified
under each combination of bias threshold (0.25, 0.5, 1.0 log₂) and minimum-tissue cutoff
(1–4 datasets), for both the full four-dataset set and the three FVB/NJ-matched datasets.

---

## Supplementary Table S4 — WASP read attrition per sample

One row per sample (all five datasets, both strains). Columns:

| Column | Description |
|--------|-------------|
| `gse` | GEO series |
| `srr` | SRR accession |
| `strain` | FVB or B6 |
| `library` | PE or SE |
| `assigned_naive` | Assigned reads in naive arm (featureCounts) |
| `assigned_wasp` | Assigned reads in WASP arm |
| `wasp_attrition_pct` | Percentage of naive assigned reads lost by WASP filtering |
| `vw2_reads` | Number of vW=2 reads in the naive BAM |
