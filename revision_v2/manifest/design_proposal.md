# Phase B-3 / B-5 — design proposal with estimates

Source of all counts: `revision_v2/manifest/manifest.tsv` (built by `revision_v2/scripts/build_manifest.py` from GEO/ENA sample attributes; `revision_v2/logs/build_manifest.log`, `B3_estimate.log`). Timing basis: old STAR logs on this machine (8 threads), median 1.30 min (naive arm) + 0.74 min (FVB arm) per million read-pairs over 18 fructose-study samples; fastp + WASP filter + featureCounts + identity QC taken as 0.6 min per million read-pairs. Disk per sample ≈ 2.4 × FASTQ size (raw + trimmed FASTQ, 3 BAMs + indices).

## Sample availability (fructose study GSE123896; animals, GEO strain attribute; 3 technical runs per animal)

| Tissue (GSE) | FVB/NJ water | FVB/NJ fructose | C57BL/6J water | C57BL/6J fructose | read length |
|---|---|---|---|---|---|
| Adipose (GSE123875) | 3 | 5 | 4 | 2 | 2 × 69 bp |
| Hypothalamus (GSE123893) | 4 | 4 | 4 | 4 | 2 × 101 bp |
| Liver (GSE123894) | 6 | 6 | 6 | 6 | 2 × 69 bp |

Sex and age are not given in the GEO/BioSample attributes of this study (columns left empty in the manifest; the original paper states males). Read length differs by tissue (hypothalamus 101 bp), not by strain, so it cannot confound the strain contrast but will be reported.

## Design 1 — main analysis (recommended): control-diet (water) animals only, model `~ strain`
- Animals: 27 (adipose 3 vs 4; hypothalamus 4 vs 4; liver 6 vs 6). 81 runs.
- Read-pairs ≈ 1,377 M; FASTQ ≈ 115.5 GB (ENA fastq_bytes).
- 3-arm wall time ≈ 61 h sequential on 8 cores (≈ 2.5 days); disk ≈ 277 GB peak if everything is kept; ≈ 150 GB if trimmed FASTQ is deleted after alignment (free space now 420 GB).
- Rationale: no treatment term; the strain contrast is interpretable directly. The old catalog rule (bias ≥ 0.5 log2 in ≥ 2 tissues) can be applied unchanged over 3 tissues.
- Caveat: adipose FVB/NJ has only 3 animals (edgeR QL F-test with 3 vs 4 is feasible but low-powered; the bias statistic itself is a within-sample arm difference and is not affected by n).

## Design 2 — sensitivity analysis: both diets, model `~ diet + strain`
- Animals: 54 (adipose 8 vs 6; hypothalamus 8 vs 8; liver 12 vs 12). 162 runs.
- Read-pairs ≈ 2,645 M; FASTQ ≈ 225 GB; wall time ≈ 116 h (≈ 5 days); disk ≈ 540 GB peak — **exceeds free space** unless processed in batches with deletion of FASTQ/BAM after counting (feasible: counts and QC tables are small).
- Use: show that the bias catalog does not depend on diet, and (reviewer 2, comment 1) demonstrate the effect of reference choice on a *within-FVB* comparison (fructose vs water in FVB/NJ, 3 tissues) — this directly answers R2-1 with real data.
- Recommendation: run Design 1 first and complete the manuscript on it; run Design 2 afterwards as an incremental batch (the 27 water animals are shared, so only the 27 fructose animals, ≈ 110 GB / ≈ 55 h, are additional) if time allows. Decision requested.

## B-5 — handling of the old datasets
- GSE135230 (microglia, mixed background) and the DBA/2J runs of GSE123875: **not used** (manifest include = N with reason).
- GSE175625 (13 B6.FVB-Tg(MMTV-PyVT) tumours, all QC "B6 pure", alt-allele fraction 0.011–0.016): two options for the human to choose.
  - **Option A (proposed): negative control.** Process the 13 tumours through the same three arms. Expectation: genome-wide bias distribution centred on 0, no catalog genes. This gives an empirical specificity check for the bias statistic that the F1 data cannot provide. Cost: 2.1 GB FASTQ, ≈ 3 h, ≈ 5 GB disk. Methods must state that these are B6-congenic animals used solely as an FVB-allele-free control; the old CNV discussion is deleted. Caveat: the residual 1–2 % FVB-homozygous sites are consistent with the congenic interval around the PyVT transgene on chromosome 6 and may produce a small local bias signal there — this is reportable, not a problem.
  - **Option B: complete exclusion.** Simpler manuscript; loses the specificity check.
- F1 (GSE200632 ×4, GSE293440 ×8): Phase A QC = F1 for all 12; reused as is (include = Y). GSE293440 tissue wording to be corrected to "cerebral (frontal) cortex" and the substrain statement removed.

## Technical-run handling
Each animal's 3 runs (same BioSample, same library) are concatenated per read (`cat R1_run1 R1_run2 R1_run3`) before fastp; this is recorded in `analysis_v2/params.yaml` in Phase C. Identity QC is run on every animal's naive BAM before the FVB arm, and the pipeline stops on any mismatch with the manifest strain.

## Decisions (2026-10-14, after Phase B STOP; author recommendation adopted)
1. **Design: both.** Primary analysis is pre-specified as `~ diet + strain` on all 54 animals (adipose 8 vs 6, hypothalamus 8 vs 8, liver 12 vs 12), which satisfies the adoption rule "biological n ≥ 4 per strain × tissue" in every tissue; the water-only subset (27 animals) is the sensitivity analysis. Both diets also allow the within-FVB diet comparison (fructose vs water, naive arm vs FVB arm) that answers reviewer 2 comment 1 with data. Execution order: water animals first, fructose animals as a second batch; tissue-wise batches; trimmed FASTQ deleted after alignment, FVB-arm and WASP-arm BAMs deleted after featureCounts and QC (naive BAM kept).
2. **PyVT: Option A (negative control).** No FVB-vs-B6 contrast exists in these 13 animals, so the bias statistic (difference of log2FC) is undefined; instead a per-sample reference-swap ratio, log2(CPM_FVBref / CPM_naive) per gene, is computed and expected to be ≈ 0 genome-wide and at catalog genes. The same metric is computed for the fructose-study B6 animals for comparison. Per-chromosome FVB-allele fractions are reported, because a congenic FVB-derived interval around the transgene insertion may give a local positive signal (internal positive control). Methods will present this as a separate metric from the main analysis.
3. **Additional dataset: verify GSE55792** (wild-type C57BL/6J and FVB/NJ cardiac polyA RNA-seq). Adopt as a fourth tissue only if sample attributes show ≥ 3 animals per strain, bulk polyA libraries, and all samples pass the sequence QC after download. Prostate / skin candidates are reserves. HMDP n = 2 sets are not used.
4. Sex and age are not recorded in GEO/BioSample for the fructose study: stated as a limitation in the manuscript.
