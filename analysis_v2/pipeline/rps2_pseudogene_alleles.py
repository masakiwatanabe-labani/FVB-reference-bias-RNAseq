#!/usr/bin/env python3
"""Phase K-1'(c): which allele (B6 or FVB/NJ) do the Rps2 processed pseudogenes ps10, ps13 and ps6 (Ensembl 110, B6 genome) carry
at the exonic FVB/NJ SNP positions of Rps2, and do they carry the FVB/NJ sequence at the four exonic indels?
Method: Rps2 mature sequence = concatenation of the merged Ensembl 110 exons (revision_v2/audit/rps2_exons.GRCm39.110.bed, + strand);
each pseudogene's genomic sequence is aligned to it with blastn (word size 11; all HSPs with e-value < 1e-50, i.e. the segments of the
processed mRNA that the pseudogene copies); at each exonic SNP inside an HSP the pseudogene base aligned to the Rps2 position is compared
with the B6 (REF) and FVB/NJ (ALT) bases; a variant outside every HSP lies in exonic sequence of the merged Ensembl exon set that is not
represented in the pseudogene ("not in pseudogene"). Indels are called the same way (FVB/NJ allele if the HSP shows the corresponding gap
or insertion at that position, B6 if it shows the reference sequence).
Output: analysis_v2/results/worked_examples/rps2_pseudogene_alleles.tsv"""
import subprocess, pandas as pd
P="/usr/local/jupyter/FVB_bias"; FA=f"{P}/ref/b6_fasta/GRCm39.primary.fa"
def fetch(reg): return "".join(subprocess.run(["samtools","faidx",FA,reg],capture_output=True,text=True,check=True).stdout.split("\n")[1:]).upper()
ex=[tuple(map(int,l.split()[1:3])) for l in open(f"{P}/revision_v2/audit/rps2_exons.GRCm39.110.bed")]
seq=""; pos2i={}
for s,e in ex:
    sq=fetch(f"17:{s+1}-{e}")
    for k in range(len(sq)): pos2i[s+1+k]=len(seq)+k
    seq+=sq
snps=[]
out=subprocess.run(["bcftools","view","-H","-r","17:24937090-24940903",f"{P}/ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf.gz"],capture_output=True,text=True,check=True).stdout
for l in out.strip().split("\n"):
    f=l.split("\t"); p=int(f[1])
    if p in pos2i: snps.append((p,f[3],f[4],f[9].split(":")[0]))
indels=[]
out=subprocess.run(["bcftools","view","-H","-r","17:24937090-24940903",f"{P}/ref/fvb_vcf/FVB_NJ.indels.GRCm39.vcf.gz"],capture_output=True,text=True,check=True).stdout
for l in out.strip().split("\n"):
    f=l.split("\t"); p=int(f[1])
    if p in pos2i and (p+len(f[3])-1) in pos2i: indels.append((p,f[3],f[4]))
gm=pd.read_csv(f"{P}/analysis_v2/results/gene_meta.tsv",sep="\t",index_col=0,dtype={"chrom":str})
import tempfile, os
tmp=tempfile.mkdtemp(); open(f"{tmp}/q.fa","w").write(">rps2\n"+seq+"\n")
rows=[]
for name in ("Rps2-ps10","Rps2-ps13","Rps2-ps6"):
    g=gm[gm.gene_name==name].iloc[0]; ps=fetch(f"{g.chrom}:{g.start}-{g.end}"); open(f"{tmp}/s.fa","w").write(">ps\n"+ps+"\n")
    out=subprocess.run(["blastn","-query",f"{tmp}/q.fa","-subject",f"{tmp}/s.fa","-word_size","11","-evalue","1e-50","-outfmt","6 qstart qend sstart send pident qseq sseq"],capture_output=True,text=True,check=True).stdout
    hsps=[l.split("\t") for l in out.strip().split("\n") if l]
    qmap={}   # mature position (0-based) -> pseudogene base ('-' if deleted); insertions in the pseudogene recorded after the preceding base
    ins={}
    for h in hsps:
        qs=int(h[0])-1; qseq,sseq=h[5],h[6]; qi=qs; last=None
        for a,b in zip(qseq,sseq):
            if a=="-": ins[last]=ins.get(last,"")+b; continue
            qmap[qi]=b; last=qi; qi+=1
    span=";".join(f"{h[0]}-{h[1]}" for h in hsps); idt=";".join(h[4] for h in hsps)
    for p,ref,alt,gt in snps:
        i=pos2i[p]; b=qmap.get(i)
        call="not in pseudogene" if b is None else "FVB/NJ" if b==alt else "B6" if b==ref else "other"
        rows.append({"pseudogene":name,"chrom":g.chrom,"start":g.start,"end":g.end,"strand":g.strand,"rps2_mature_segments":span,"identity_pct":idt,
                     "variant":f"SNP 17:{p} {ref}>{alt}","FVB_NJ_GT":gt,"rps2_mature_pos":i+1,"pseudogene_base":b or "","allele_carried":call})
    for p,ref,alt in indels:
        i=pos2i[p]; inside=all((i+k) in qmap for k in range(len(ref)))
        if not inside: call="not in pseudogene"; b=""
        else:
            if len(alt)>len(ref): b=ins.get(i+len(ref)-1,""); call="FVB/NJ" if b==alt[len(ref):] else "B6" if b=="" else "other"
            else: dels=[qmap[i+k] for k in range(1,len(ref))]; b="".join(dels); call="FVB/NJ" if set(dels)=={"-"} else "B6" if "-" not in dels else "other"
        rows.append({"pseudogene":name,"chrom":g.chrom,"start":g.start,"end":g.end,"strand":g.strand,"rps2_mature_segments":span,"identity_pct":idt,
                     "variant":f"indel 17:{p} {ref}>{alt}","FVB_NJ_GT":"1/1","rps2_mature_pos":i+1,"pseudogene_base":b,"allele_carried":call})
df=pd.DataFrame(rows); df.to_csv(f"{P}/analysis_v2/results/worked_examples/rps2_pseudogene_alleles.tsv",sep="\t",index=False)
print(f"Rps2 mature length {len(seq)}; exonic SNPs {len(snps)}; exonic indels {len(indels)}")
pd.set_option("display.width",220); print(df[["pseudogene","rps2_mature_segments","variant","rps2_mature_pos","pseudogene_base","allele_carried"]].to_string(index=False))
print(df.groupby(["pseudogene","allele_carried"]).size().unstack(fill_value=0))
