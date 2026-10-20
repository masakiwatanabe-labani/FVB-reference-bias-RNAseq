# GSE55792 (heart) — sequence-based identity check of the "C57BL/6J" group (2026-10-14)

Sample heart_B6_C57_M_12wk_mRNA_3 (GSM1346089, GSE55788 "C57BL/6J wild type, 12 wk male, whole heart", SRR1187093, SE 50 bp) was processed with `analysis_v2/pipeline/run_sample.sh` (download → fastp → STAR naive → ASEReadCounter at 203,510 FVB/NJ hom-alt exonic sites).

Result (`analysis_v2/data/heart/heart_B6_C57_M_12wk_mRNA_3/qc/identity.tsv`): 3,693 sites with depth ≥ 20; mean alt fraction 0.187; frac_hi 0.162; frac_mid 0.052; frac_lo 0.786 → verdict **needs review**, expected **B6 pure** → pipeline stopped (rule 0.4-4).

Genomic distribution of the homozygous-alt sites (20-Mb bins, `ase.tsv`): concentrated in large blocks — chr2 100–140 Mb (69–95 % hom-alt), chr3 20–40 Mb, chr8 40–100 Mb (77–100 %), chr10 20–40 Mb, chr11 40–100 Mb (61–80 %), chr12 0–20 Mb, chr17 60–80 Mb — while chr1, 5, 7, 13 are ≥ 99 % reference. This is the pattern of an animal carrying large homozygous non-B6 (FVB-allele-sharing) segments, i.e. an incompletely backcrossed / mixed background, not C57BL/6J. The calibration B6 controls (fructose study) are ≥ 99.8 % reference at the same sites.

Consequence: this animal cannot be used as a B6 comparator. The remaining three GSE55788 "C57BL/6J" animals are being checked the same way (`analysis_v2/logs/batch_chain.log`). If they fail too, GSE55792 is dropped from the manifest (include=N, reason "sequence QC: not B6") and the heart tissue is not added; the prostate / skin candidates remain as reserves (human decision). The FVB/NJ heart animals will not be downloaded unless a B6 group is available.

## Update 2026-10-14 12:25 — all four GSE55788 "C57BL/6J" animals fail

| sample | GSM | n_sites | mean alt fraction | frac_hi | frac_mid | frac_lo | verdict |
|---|---|---|---|---|---|---|---|
| heart_B6_C57_M_12wk_mRNA_1 | GSM1346087 | 5,124 | 0.205 | 0.180 | 0.055 | 0.765 | needs review |
| heart_B6_C57_M_12wk_mRNA_2 | GSM1346088 | 5,038 | 0.196 | 0.149 | 0.103 | 0.748 | needs review |
| heart_B6_C57_M_12wk_mRNA_3 | GSM1346089 | 3,693 | 0.187 | 0.162 | 0.052 | 0.786 | needs review |
| heart_B6_C57_M_12wk_mRNA_4 | GSM1346090 | 5,504 | 0.207 | 0.180 | 0.059 | 0.761 | needs review |

The homozygous-alt sites sit in the same large chromosomal blocks in all four animals (see run log below), i.e. the colony carries fixed non-B6 segments. **Decision (2026-10-14): GSE55792 is not used** (manifest include=N for all 111 samples with this reason); the FVB/NJ heart samples were not downloaded. The four naive BAMs and QC tables are kept under `analysis_v2/data/heart/` as evidence (trimmed FASTQ deleted). The strain-comparison analysis proceeds with the three fructose-study tissues (adipose, hypothalamus, liver), as the instructions allow when no additional dataset qualifies (B-4). Reserve candidates (GSE86315 skin controls; prostate SRP303740/SRP430873) are documented in `candidate_datasets.md`; GSE86315 attributes are checked below for a possible later batch.
