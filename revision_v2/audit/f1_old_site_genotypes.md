# H-1a — Which sites did the old F1 ASE analysis evaluate? (2026-10-20)

**Question.** The Results of v2 reported, for the same 12 F1 animals, a mean FVB-allele fraction of 0.474–0.481 in the sample-identity QC (203,510 FVB/NJ-homozygous exonic sites, depth ≥ 20) but 0.333 in the Fig. 2 ASE analysis (274–872 sites per animal). Fig. 2 re-used the per-site GATK ASEReadCounter output of the original (v1) run.

**Inputs.** The 12 old naive per-site tables (`analysis_v2/results/f1_ase/GSE*/<SRR>.naive.ase.tsv`, byte-identical copies of the legacy output; md5 in `analysis_v2/results/f1_ase/copied_from_legacy.md5`); the VCF the old run passed to ASEReadCounter — the legacy logs name no VCF, but the old pipeline's STAR/WASP step used `ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf` (bcftools `-s FVB_NJ -f PASS --min-ac 1`, **no genotype filter**, 5,960,448 records), and every evaluated site is in that file (below); the QC site set `analysis_v2/qc/fvbnj_homalt_exonic.hetGT.vcf.gz` (203,510 sites).

**Method.** Union of the (contig, position) pairs of the 12 old naive tables → `revision_v2/audit/f1_old/old_sites.bed`; `bcftools query -R` against the two VCFs → `revision_v2/audit/f1_old/old_sites.snps_vcf_GT.tsv`.

## Result

| item | value |
|---|---|
| distinct sites evaluated by the old run (union over 12 naive tables) | 94,351 |
| of which FVB_NJ genotype **0/1** in `FVB_NJ.snps.GRCm39.vcf.gz` | **94,351 (100%)** |
| of which genotype 1/1 or other | 0 |
| overlap with the 203,510 QC sites (FVB_NJ 1/1, PASS, biallelic, autosomal, exonic) | **0** |
| sites per sample evaluated (raw) / passing the old filter (total ≥ 20, ref ≥ 3, alt ≥ 3) | 10,490–29,218 / 274–872 (legacy `f1_ase_combined/f1_ase_all_samples.tsv`) |
| fraction of raw evaluated sites with altCount = 0 (legacy summary `alt0_pct`) | 73–78% |

**Interpretation.** GATK ASEReadCounter evaluates only sites whose VCF genotype is heterozygous. The old run supplied the unfiltered FVB/NJ SNP VCF, in which the sites that matter for an F1 — those where FVB/NJ is homozygous for the alternative allele — carry genotype 1/1 and were therefore **skipped**; the sites that were evaluated are exclusively the 94,351 positions at which the Mouse Genomes Project calls FVB_NJ **heterozygous**, i.e. positions where the FVB/NJ inbred genotype is uncertain (low-quality or copy-number-affected calls). At such sites an F1 animal is not expected to show a 0.50 allele fraction; most of them showed no alternative allele at all (73–78% altCount = 0), and the ≥ 3/≥ 3 filter kept only the minority with some alternative reads, whose mean (0.317–0.348 per sample) does not measure reference bias at genuine heterozygous sites. The v2 identity QC, which set the genotype of the 203,510 homozygous-alt sites to 0/1 so that every site is evaluated, measured 0.474–0.481 on the same BAMs.

**Consequence.** The old F1 values (0.333 naive, 0.386 WASP, "17 percentage points") are an artefact of site selection, not a measurement of reference bias. H-1b recomputes the ASE at the QC site set. Because the old naive/WASP BAMs of the F1 animals were deleted with the legacy data (`revision_v2/audit/legacy_data_inventory_before_deletion.tsv`; no copy under `analysis_v2/data`), the 12 animals are re-downloaded and re-aligned with the v2 pipeline (`analysis_v2/pipeline/run_sample.sh`, naive arm; WASP BAM regenerated with `samtools view -e '![vW] || [vW]==1'`) before the recount — see `analysis_v2/pipeline/f1_ase_rerun.sh` and CHANGELOG.
