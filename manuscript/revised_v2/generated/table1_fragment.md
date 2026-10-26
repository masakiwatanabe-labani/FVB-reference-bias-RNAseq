# Table 1. RNA-seq datasets used in this study


## Three-arm strain comparison datasets (FVB/NJ versus C57BL/6J)

| GEO accession | Tissue | Cohort | *n* FVB/NJ (conditions) | *n* C57BL/6J (conditions) | Library | Read length | Runs per sample | Sequence identity QC passed |
|---|---|---|---|---|---|---|---|---|
| GSE123875 | Mesenteric adipose | Fructose study (PRJNA510125; SubSeries of GSE123896) | 8 (fructose ×5; water ×3) | 6 (fructose ×2; water ×4) | Paired-end | 69 bp | 3 | 14/14 |
| GSE123893 | Hypothalamus | Fructose study (PRJNA510129; SubSeries of GSE123896) | 8 (fructose ×4; water ×4) | 8 (fructose ×4; water ×4) | Paired-end | 101 bp | 3 | 16/16 |
| GSE123894 | Liver | Fructose study (PRJNA510133; SubSeries of GSE123896) | 12 (fructose ×6; water ×6) | 12 (fructose ×6; water ×6) | Paired-end | 69 bp | 3 | 24/24 |
| GSE86315 | Back skin | Imiquimod strain panel (PRJNA341396) | 8 (IMQ F ×2; IMQ M ×2; control F ×2; control M ×2) | 8 (IMQ F ×2; IMQ M ×2; control F ×2; control M ×2) | Single-end | 52 bp | 1 | 16/16 |

## Negative-control dataset (no FVB/NJ-versus-C57BL/6J contrast)

| GEO accession | Tissue | Strain (GEO attribute) | *n* | Conditions | Library | Read length | Sequence identity QC |
|---|---|---|---|---|---|---|---|
| GSE175625 | Mammary tumour (MMTV-PyVT) | B6.FVB-Tg(MMTV-PyVT)634Mul/LellJ (C57BL/6 congenic) | 13 | Methylparaben ×4; Propylparaben ×6; Vehicle ×3 | Single-end | 101 bp | 13/13 (expected: C57BL/6) |

## F1 hybrid datasets used for allele-specific expression analysis (Fig. 2)

| GEO accession | Tissue | Cross (GEO attribute) | *n* (F1 samples) | Library | Read length | Sequence identity QC |
|---|---|---|---|---|---|---|
| GSE293440 | Cerebral (frontal) cortex | B6xFVB | 8 | Paired-end | 150 bp | 8/8 (expected: F1) |
| GSE200632 | Primary keratinocytes (LSL-KrasG12D, Cre-transduced) | C57BL/6:FVB/N (F1:50/50) | 4 | Paired-end | 151 bp | 4/4 (expected: F1) |

Adipose, hypothalamus and liver samples of the fructose study were collected from the same animals; *n* per tissue counts samples (fructose study: FVB/NJ 28 samples from 15 animals, C57BL/6J 26 samples from 15 animals; all strain-comparison samples: FVB/NJ 36 samples from 23 animals, C57BL/6J 34 samples from 23 animals; Supplementary Table S1). The PyVT tumours are one primary tumour per mouse; GSE200632 *n* counts keratinocyte cultures from newborn pups (number of donor animals not recorded).
Sequence identity QC: allele fractions at FVB/NJ-homozygous exonic SNP sites (203,510; a fixed 10% subset for all but the F1 samples) in the B6-reference alignment (see *Sample identity verification*); full per-sample results in Supplementary Table S1.
Sex and age are not recorded in the GEO/BioSample attributes of the fructose study; the imiquimod study records sex (M/F). All FVB/NJ animals are from the Jackson Laboratory substrain, matched to the FVB_NJ variant set.
