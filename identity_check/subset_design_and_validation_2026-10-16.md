# Speed-up Step 3 — strain-identity QC on a fixed 10 % site subset (2026-10-16)

## What was changed

| | before | after (from hypothalamus_B6_CWHyp4 onward) |
|---|---|---|
| sites | `analysis_v2/qc/fvbnj_homalt_exonic.hetGT.vcf.gz`, 203,510 FVB/NJ hom-alt, PASS, biallelic, autosomal, exonic SNPs (MGP REL-2021, GT rewritten to 0/1) | `fvbnj_homalt_exonic.hetGT.sub10pct.seed20261016.vcf.gz`: 20,351 sites = `random.Random(20261016).sample(range(203510), 20351)` of the same file, order preserved |
| GATK traversal | whole BAM | `-L fvbnj_homalt_exonic.hetGT.sub10pct.seed20261016.intervals.bed` (the subset sites merged within 1 kb: 13,507 intervals, 2.11 Mb). `-L` only limits where the walker goes; the per-site pileup and counts are unchanged |
| ASEReadCounter options | `--min-depth-of-non-filtered-base 20 --min-mapping-quality 20 --min-base-quality 20` | unchanged |
| metrics and verdict rules (`identity_verdict.py`) | mean_af, frac_hi (>0.9), frac_mid, frac_lo (<0.1) over sites with depth ≥ 20; rules fixed before Phase A | unchanged |

md5: `analysis_v2/qc/qc_sites_subsample.md5` (subset VCF 6d4b3b7e…, full VCF 9e650aea…, intervals BED).

## Validation 1 — all 13 samples evaluated so far, full set vs subset (subset = full-set `ase.tsv` filtered to the subset sites)

| sample | n sites full | mean_af full | frac hi/mid/lo full | verdict full | n sites 10 % | mean_af 10 % | frac hi/mid/lo 10 % | verdict 10 % | same |
|---|---|---|---|---|---|---|---|---|---|
| adipose_B6_CWM5 | 25225 | 0.007 | 0.001/0.006/0.993 | B6 pure | 2590 | 0.007 | 0.002/0.004/0.994 | B6 pure | Y |
| adipose_B6_CWM6 | 27413 | 0.004 | 0.001/0.002/0.997 | B6 pure | 2765 | 0.003 | 0.001/0.001/0.998 | B6 pure | Y |
| adipose_B6_CWM7 | 18549 | 0.002 | 0.001/0.001/0.999 | B6 pure | 1868 | 0.003 | 0.002/0.000/0.998 | B6 pure | Y |
| adipose_B6_CWM8 | 35058 | 0.003 | 0.002/0.001/0.997 | B6 pure | 3467 | 0.003 | 0.002/0.001/0.997 | B6 pure | Y |
| adipose_FVB_FWM1 | 20056 | 0.991 | 0.989/0.004/0.007 | FVB/NJ pure | 2049 | 0.990 | 0.988/0.004/0.008 | FVB/NJ pure | Y |
| adipose_FVB_FWM7 | 17515 | 0.992 | 0.990/0.003/0.007 | FVB/NJ pure | 1760 | 0.991 | 0.989/0.004/0.007 | FVB/NJ pure | Y |
| adipose_FVB_FWM9 | 30685 | 0.991 | 0.989/0.005/0.005 | FVB/NJ pure | 3103 | 0.992 | 0.991/0.004/0.005 | FVB/NJ pure | Y |
| heart_B6_C57_M_12wk_mRNA_1 (GSE55788) | 5124 | 0.205 | 0.180/0.055/0.765 | needs review | 519 | 0.203 | 0.175/0.062/0.763 | needs review | Y |
| heart_B6_C57_M_12wk_mRNA_2 | 5038 | 0.196 | 0.149/0.103/0.748 | needs review | 518 | 0.187 | 0.141/0.100/0.759 | needs review | Y |
| heart_B6_C57_M_12wk_mRNA_3 | 3693 | 0.187 | 0.162/0.052/0.786 | needs review | 366 | 0.181 | 0.158/0.052/0.790 | needs review | Y |
| heart_B6_C57_M_12wk_mRNA_4 | 5504 | 0.207 | 0.180/0.059/0.761 | needs review | 553 | 0.196 | 0.166/0.067/0.767 | needs review | Y |
| hypothalamus_B6_CWHyp1 | 34398 | 0.002 | 0.001/0.001/0.998 | B6 pure | 3482 | 0.002 | 0.001/0.001/0.998 | B6 pure | Y |
| hypothalamus_B6_CWHyp2 | 34550 | 0.002 | 0.001/0.001/0.999 | B6 pure | 3518 | 0.002 | 0.001/0.002/0.997 | B6 pure | Y |

Largest differences: mean_af 0.010 (heart 2), frac_hi 0.014 (heart 4), frac_mid 0.008; all 13 verdicts identical, including the four "needs review" heart samples whose metrics sit nearest to a rule boundary. 2,000–3,500 covered sites per animal remain (tissue-dependent).

The subset-based tables were written as `qc/identity.sub10pct.tsv` (and `qc/<sample>.ase.sub10pct.tsv`) next to the untouched full-set `identity.tsv` / `ase.tsv` of these 13 samples, so that the Supplementary "Sample identity verification" table can report one site set for every animal; the full-set values stay available as evidence (this file, CHANGELOG).

## Validation 2 — real run with the subset VCF + intervals (adipose_B6_CWM5, 48 M pairs)

| run | sites | `-L` | wall time | output |
|---|---|---|---|---|
| production (2026-10-14) | 203,510 | none | 30 min | `qc/adipose_B6_CWM5.ase.tsv` (25,225 sites with depth ≥ 20) |
| subset + merged intervals (new) | 20,351 | 13,507 intervals | **8.1 min** | **byte-identical** to the production output filtered to the subset sites (2,590 sites) |
| full set + sites merged within 10 kb (tried, not adopted) | 203,510 | 20,191 intervals | ≈ 13 min (see note) | — |
| full set + one interval per site (tried, not adopted) | 203,510 | 203,510 | slower than no `-L` (≈ 70 sites/s; stopped) | — |

ASEReadCounter time is dominated by walking covered loci, not by the number of sites, which is why the subset alone gains little and the interval restriction is what makes it fast; a single-base interval per site is slower because of per-interval index seeks. The subset + 1-kb-merged intervals is 3.7× faster than production, and with Step 1 the step overlaps with STAR anyway.

## Where it is recorded

- `analysis_v2/params.yaml`: `qc_sites_full`, `qc_sites`, `qc_sites_intervals`.
- `run_sample.sh`: `SITES`, `SITES_BED`, `-L $SITES_BED` (step 5). Old copy: `revision_v2/speedup/run_sample.sh.before_step1`.
- Methods (Phase D), Sample identity verification: "Allele fractions were evaluated with GATK ASEReadCounter (minimum depth 20, MAPQ ≥ 20, base quality ≥ 20) at a fixed random subset of 20,351 of the 203,510 exonic FVB/NJ homozygous-alternative SNPs (seed 20261016); the first ten animals and the four GSE55788 heart samples were additionally evaluated at the full set with identical verdicts (Supplementary Table …)."
