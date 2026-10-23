#!/usr/bin/env python3
"""Fixed 10% subset of the sample-identity check sites (as used from 2026-10-16 for every animal except the F1 hybrids).
Full set: the 203,510 FVB/NJ homozygous-alternative, PASS, biallelic, autosomal, exonic SNPs of MGP REL-2021
(revision_v2/scripts/A2_build_qc_sites.sh -> fvbnj_homalt_exonic.vcf.gz), with the genotype rewritten to 0/1 for GATK
ASEReadCounter:  bcftools +setGT -Oz -o fvbnj_homalt_exonic.hetGT.vcf.gz -- fvbnj_homalt_exonic.vcf.gz -t a -n c:0/1
Subset: random.Random(20261016).sample(range(203510), 20351), records kept in file order.
Traversal intervals for GATK -L (limit where the walker goes; per-site counts unchanged): the subset sites merged with
  bedtools merge -d 999  (reproduces identity_intervals_sub10pct.seed20261016.bed byte for byte).
Usage: make_identity_subset.py fvbnj_homalt_exonic.hetGT.vcf.gz OUT_PREFIX
Writes OUT_PREFIX.vcf (subset VCF, header + records), OUT_PREFIX.sites.tsv (CHROM POS REF ALT)."""
import gzip, random, sys
full_vcf, out = sys.argv[1], sys.argv[2]
hdr, rec = [], []
for l in gzip.open(full_vcf, "rt"):
    (hdr if l.startswith("#") else rec).append(l)
assert len(rec) == 203510, len(rec)
keep = sorted(random.Random(20261016).sample(range(len(rec)), 20351))
with open(f"{out}.vcf", "w") as v, open(f"{out}.sites.tsv", "w") as s:
    v.writelines(hdr); s.write("CHROM\tPOS\tREF\tALT\n")
    for i in keep:
        v.write(rec[i]); f = rec[i].split("\t", 5); s.write(f"{f[0]}\t{f[1]}\t{f[3]}\t{f[4]}\n")
print(f"{len(keep)} sites -> {out}.vcf, {out}.sites.tsv")
