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

## (b) Reads at Rps2 and its processed pseudogenes in skin (GSE86315), per strain — revised 2026-10-22 (Phase K-1′(a))

**Definition used everywhere in the article** (Results, Supplementary Note 2, numbers.tsv keys `rps2_skin_*`, `rps2_pseudogenes_skin_*`, `rps2mm_*`): shared-denominator mean CPM — each arm and counting mode divided by the naive-arm unique-only library size of the same animal (the normalisation of Fig. 6B), mean over the 8 FVB/NJ and 8 C57BL/6J skin animals — for Rps2 and the **three** processed pseudogenes that carry the multi-mapping signal, Rps2-ps10, Rps2-ps13 and Rps2-ps6. Source: `analysis_v2/results/worked_examples/rps2_multimap_cpm.tsv` (`analysis_v2/pipeline/rps2_multimap_cpm.py`).

| gene | reference | FVB/NJ unique | FVB/NJ −M −−fraction | C57BL/6J unique | C57BL/6J −M −−fraction |
|---|---|---|---|---|---|
| Rps2 | B6 | 3.71 | 13.00 | 209.76 | 232.82 |
| Rps2 | FVB/NJ | 3.67 | 41.84 | 5.06 | 19.90 |
| ps10 / ps13 / ps6 | B6 | 0.05 / 0.13 / 0.10 | 49.42 / 44.41 / 7.64 | 0.01 / 0.08 / 0.11 | 9.68 / 7.24 / 2.49 |
| ps10 / ps13 / ps6 | FVB/NJ | 0.05 / 0.13 / 0.09 | 37.83 / 32.58 / 6.23 | 0.01 / 0.11 / 0.12 | 14.69 / 9.48 / 6.43 |

Three pseudogenes summed, −M −−fraction: FVB/NJ 101.5 (B6 reference) and 76.6 (FVB/NJ reference); C57BL/6J 19.4 and 30.6. Rps2 + three pseudogenes, FVB/NJ: 114.5 and 118.5.

**Cause of the earlier numbers.** The values quoted in the Phase H-2b/J text (111 versus 29 under the B6 reference; 87 versus 42 under the FVB/NJ reference) were the sums over **all 12** annotated Rps2 processed pseudogenes (Rps2-ps2 … ps13; table `rps2_pseudogene_counts_skin.tsv`, kept as the record), whereas the sentence named three. Denominator and averaging were the same; the difference (≈ 9–12 CPM) is the fractional counts on the other nine pseudogenes (mostly Rps2-ps8 and -ps7). The article now uses the three-pseudogene values throughout.

Rps2 itself rises under fractional counting in both strains and on both references, by a similar log₂ amount (FVB/NJ/C57BL/6J log₂ ratio −5.82 → −4.16 on the B6 reference, −0.46 → +1.07 on the FVB/NJ reference), so the bias (difference of the two) is 5.36 unique-only and 5.23 fractional — the values of Supplementary Table S2 — and the residual fraction stays 0.98.

## (c) Alleles carried by the pseudogenes (Phase K-1′(c); `analysis_v2/results/worked_examples/rps2_pseudogene_alleles.tsv`)

blastn of each pseudogene against the Rps2 sequence (merged Ensembl 110 exons): the pseudogenes copy the processed mRNA segments 365/368–742, 883–1224 and 1396–1575/1616 (99.1–99.5% identity). Only 5 of the 13 exonic FVB/NJ SNPs of Rps2 lie in these segments (the other 8, and all 4 exonic indels, lie in exon sequence of other Ensembl transcripts that the pseudogenes do not contain). At these 5 SNPs Rps2-ps10 and Rps2-ps13 carry the FVB/NJ allele at all 5 and Rps2-ps6 at 3. On the SNP-patched FVB/NJ reference, Rps2, ps10 and ps13 are therefore identical at every variant position within the shared segments, which is consistent with the near-equal division of FVB/NJ reads among them under fractional counting; on the B6 reference, FVB/NJ reads match ps10 and ps13 better than Rps2.
