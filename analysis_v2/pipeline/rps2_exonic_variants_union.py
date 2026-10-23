#!/usr/bin/env python3
"""Phase M-3 D2: Rps2 exonic FVB/NJ variants with a stated, reproducible exon definition.
Exonic = within the union of the exons of all Ensembl 110 Rps2 transcripts (the featureCounts meta-feature of gene-level counting),
taken directly from ref/b6_gtf/GRCm39.110.gtf (merged). Variants: FVB_NJ homozygous-alternative SNPs (FVB_NJ.snps.GRCm39.vcf.gz) and
indels (FVB_NJ.indels.GRCm39.vcf.gz) overlapping that union.
Copied part: the blastn HSP spans already used for rps2_pseudogene_alleles.tsv (column rps2_mature_segments, positions in the mature
sequence built from revision_v2/audit/rps2_exons.GRCm39.110.bed), mapped back to genomic positions. A variant is inside the copied part of a
pseudogene if all its reference bases lie inside that pseudogene's HSP spans. Cross-check: each pseudogene is also aligned with blastn
(same word size and e-value) to the Rps2 gene body itself, which includes the retained-intron and other exon segments of the union; any
union variant inside a genomic HSP is reported in column genomic_hsp_check.
For SNPs inside the copied part the pseudogene allele is taken from rps2_pseudogene_alleles.tsv; a SNP inside the copied part that is not
in that file would be called here from the genomic HSP (none expected: the HSP spans lie within the old exon set).
Output: analysis_v2/results/worked_examples/rps2_exonic_variants_union.tsv. Usage: rps2_exonic_variants_union.py"""
import subprocess, tempfile, re
import pandas as pd
P = "/usr/local/jupyter/FVB_bias"; FA = f"{P}/ref/b6_fasta/GRCm39.primary.fa"; GTF = f"{P}/ref/b6_gtf/GRCm39.110.gtf"
OUT = f"{P}/analysis_v2/results/worked_examples/rps2_exonic_variants_union.tsv"
def fetch(reg): return "".join(subprocess.run(["samtools", "faidx", FA, reg], capture_output=True, text=True, check=True).stdout.split("\n")[1:]).upper()
# union of Rps2 exons (all transcripts)
ex = []
for l in subprocess.run(["grep", "-P", 'gene_name "Rps2";', GTF], capture_output=True, text=True, check=True).stdout.splitlines():
    f = l.split("\t")
    if f[2] == "exon": ex.append((int(f[3]) - 1, int(f[4])))
    if f[2] == "gene": gchrom, gs, ge = f[0], int(f[3]), int(f[4])
ex.sort(); un = []
for s, e in ex:
    if un and s <= un[-1][1]: un[-1] = (un[-1][0], max(un[-1][1], e))
    else: un.append((s, e))
def in_union(p0, p1): return any(p1 > s and p0 < e for s, e in un)       # half-open [p0,p1)
def vcf(path):
    out = subprocess.run(["bcftools", "view", "-H", "-i", 'GT="AA"', "-r", f"{gchrom}:{gs}-{ge}", path], capture_output=True, text=True, check=True).stdout
    return [l.split("\t") for l in out.strip().split("\n") if l]
var = []
for f in vcf(f"{P}/ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf.gz"):
    p = int(f[1])
    if in_union(p - 1, p): var.append(dict(type="SNP", pos=p, ref=f[3], alt=f[4], length=1, ref_span=(p - 1, p)))
for f in vcf(f"{P}/ref/fvb_vcf/FVB_NJ.indels.GRCm39.vcf.gz"):
    p = int(f[1]); r, a = f[3], f[4]
    if in_union(p - 1, p - 1 + len(r)): var.append(dict(type="insertion" if len(a) > len(r) else "deletion", pos=p, ref=r, alt=a, length=abs(len(a) - len(r)), ref_span=(p - 1, p - 1 + len(r))))
# old mature coordinates -> genomic, and the HSP spans already used
old = [tuple(map(int, l.split()[1:3])) for l in open(f"{P}/revision_v2/audit/rps2_exons.GRCm39.110.bed")]
i2pos = []
for s, e in old: i2pos += list(range(s + 1, e + 1))   # mature index (0-based) -> 1-based genomic position
al = pd.read_csv(f"{P}/analysis_v2/results/worked_examples/rps2_pseudogene_alleles.tsv", sep="\t")
PSG = ("Rps2-ps10", "Rps2-ps13", "Rps2-ps6"); copied = {}
for ps in PSG:
    seg = al[al.pseudogene == ps].rps2_mature_segments.iloc[0]; pos = set()
    for a_, b_ in (map(int, x.split("-")) for x in seg.split(";")): pos |= {i2pos[k - 1] for k in range(a_, b_ + 1)}
    copied[ps] = pos
# genomic cross-check: pseudogene vs the Rps2 gene body
gm = pd.read_csv(f"{P}/analysis_v2/results/gene_meta.tsv", sep="\t", index_col=0, dtype={"chrom": str})
tmp = tempfile.mkdtemp(); body = fetch(f"{gchrom}:{gs}-{ge}"); open(f"{tmp}/q.fa", "w").write(">rps2\n" + body + "\n")
ghsp = {}
for ps in PSG:
    g = gm[gm.gene_name == ps].iloc[0]; open(f"{tmp}/s.fa", "w").write(">ps\n" + fetch(f"{g.chrom}:{g.start}-{g.end}") + "\n")
    out = subprocess.run(["blastn", "-query", f"{tmp}/q.fa", "-subject", f"{tmp}/s.fa", "-word_size", "11", "-evalue", "1e-50", "-outfmt", "6 qstart qend"], capture_output=True, text=True, check=True).stdout
    ghsp[ps] = [(gs + int(a) - 1, gs + int(b) - 1) for a, b in (l.split("\t") for l in out.strip().split("\n") if l)]
rows = []
for v in var:
    span = range(v["ref_span"][0] + 1, v["ref_span"][1] + 1)
    row = dict(type=v["type"], chrom=gchrom, pos=v["pos"], ref=v["ref"], alt=v["alt"], length_bp=v["length"],
               in_old_exon_bed="Y" if all(any(s < q <= e for s, e in old) for q in span) else "N")
    for ps in PSG:
        inside = all(q in copied[ps] for q in span); row[f"in_copied_part_{ps}"] = "Y" if inside else "N"
        if v["type"] == "SNP":
            hit = al[(al.pseudogene == ps) & (al.variant == f"SNP 17:{v['pos']} {v['ref']}>{v['alt']}")]
            row[f"allele_{ps}"] = (hit.allele_carried.iloc[0] if len(hit) else "NEW SNP INSIDE COPIED PART — call needed") if inside else "not in pseudogene"
        else:
            hit = al[(al.pseudogene == ps) & (al.variant == f"indel 17:{v['pos']} {v['ref']}>{v['alt']}")]
            row[f"allele_{ps}"] = (hit.allele_carried.iloc[0] if len(hit) else "NEW INDEL INSIDE COPIED PART — call needed") if inside else "not in pseudogene"
        row[f"genomic_hsp_{ps}"] = "Y" if any(a <= q <= b for a, b in ghsp[ps] for q in span) else "N"
    row["in_copied_part_any"] = "Y" if any(row[f"in_copied_part_{ps}"] == "Y" for ps in PSG) else "N"
    row["genomic_hsp_check"] = "consistent" if all((row[f"genomic_hsp_{ps}"] == "Y") == (row[f"in_copied_part_{ps}"] == "Y") for ps in PSG) else "DIFFERS"
    rows.append(row)
df = pd.DataFrame(rows); df.to_csv(OUT, sep="\t", index=False)
print("union exons:", " ".join(f"{gchrom}:{s+1}-{e}" for s, e in un)); print("genomic HSPs:", {k: v for k, v in ghsp.items()})
pd.set_option("display.width", 250); print(df.to_string(index=False))
print("exonic SNPs", int((df.type == "SNP").sum()), "| exonic indels", int((df.type != "SNP").sum()),
      "| SNPs inside copied part", int(((df.type == "SNP") & (df.in_copied_part_any == "Y")).sum()), "| indels inside copied part", int(((df.type != "SNP") & (df.in_copied_part_any == "Y")).sum()))
