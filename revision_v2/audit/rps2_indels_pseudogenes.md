# Rps2 fact check (Phase H-2b follow-up, user item 4; 2026-10-21)

## (a) FVB/NJ indels in Rps2 (ENSMUSG00000044533; chr17:24,937,090–24,940,903, Ensembl 110)
Source: ref/fvb_vcf/FVB_NJ.indels.GRCm39.vcf.gz (MGP REL-2021 indels, FVB_NJ sample; FILTER '.' for every record, hence omitted by g2gtools --pass). Exons: revision_v2/audit/rps2_exons.GRCm39.110.bed (6 merged exons). Command: bcftools view -r 17:24937090-24940903 | bedtools intersect -u.

Indels in the gene body (all genotype 1/1 or 1|1):
- 17:24938320 T>TTAAA GT=1|1
- 17:24938637 CA>C GT=1|1
- 17:24939228 CTCT>C GT=1/1
- 17:24940293 TA>T GT=1|1
- 17:24940307 CTT>C GT=1|1
- 17:24940352 T>TTTTG GT=1/1
- 17:24940378 CTG>C GT=1/1
- 17:24940417 CAG>C GT=1/1
- 17:24940464 G>GC GT=1|1
- 17:24940465 T>TGGGATCAAAGGC GT=1|1
- 17:24940562 C>CTT GT=1/1
- 17:24940768 T>TC GT=1/1

Within exons (4):
- 17:24939228 CTCT>C
- 17:24940464 G>GC
- 17:24940465 T>TGGGATCAAAGGC
- 17:24940562 C>CTT

Counts: gene body 12 indels (homozygous alt), exonic 4 (including a 12-bp insertion at 24,940,465 and a 3-bp deletion at 24,939,228); SNPs: 53 in the gene body (FVB_NJ.snps.GRCm39.vcf.gz), 13 within exons. None of the indels is in the SNP-only personalised reference.

## (b) Reads at Rps2 and its 12 annotated processed pseudogenes (Rps2-ps2…ps13) in skin (GSE86315), per strain
Source: analysis_v2/results/counts/skin.{naive,wasp,fvb}.{unique,multi}.tsv (featureCounts unique-only and -M --fraction), CPM on the naive-arm library size (shared denominator), mean over 8 FVB/NJ and 8 C57BL/6J animals; table revision_v2/audit/rps2_pseudogene_counts_skin.tsv.

| arm | Rps2 FVB/NJ | Rps2 C57BL/6J | 12 pseudogenes summed FVB/NJ | summed C57BL/6J |
|---|---|---|---|---|
| naive.unique | 3.7 | 209.8 | 1.3 | 1.2 |
| naive.multi | 13.0 | 232.8 | 110.7 | 28.7 |
| wasp.unique | 3.5 | 4.8 | 1.3 | 1.2 |
| fvb.unique | 3.7 | 5.1 | 1.3 | 1.3 |
| fvb.multi | 41.8 | 19.9 | 87.3 | 42.2 |

Unique-only assignment to the pseudogenes is negligible in both strains (≈ 1.2–1.3 CPM in total, every arm). With multi-mapping reads included, FVB/NJ animals carry 111 CPM on the pseudogenes under the B6 reference (Rps2-ps10 49, Rps2-ps13 44, Rps2-ps6 8; C57BL/6J 29) and 13 CPM on Rps2 itself (C57BL/6J 233); under the FVB/NJ reference 87 CPM on the pseudogenes (C57BL/6J 42) and 42 on Rps2 (C57BL/6J 20). Reading: the FVB/NJ Rps2-derived reads are present in the alignments but are multi-mappers shared between Rps2 and processed pseudogenes under both references, so unique-only counting assigns them to neither; they are not simply unaligned.
