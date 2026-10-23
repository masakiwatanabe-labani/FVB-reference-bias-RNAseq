#!/usr/bin/env python3
"""Supplementary figure: genome-wide distribution of homozygous FVB/NJ-allele sites in the four GSE55788 "wild-type C57BL/6J"
heart samples, with a verified C57BL/6J and a verified FVB/NJ animal of the fructose study as references (Phase M-3 D1b: the animal of
each strain with the most sites of depth >= 20 among the v2-pipeline animals evaluated at the full site set; candidates.tsv).
Heatmap: rows = samples, columns = 20-Mb genomic bins (autosomes 1-19), colour = fraction of sites (depth >= 20) with
alt-allele fraction > 0.9 (sequential single hue; bins with < 10 sites are blank). Vector PDF + 600-dpi TIFF.
Inputs: analysis_v2/results/identity/full_site_tables/ only (per-site tables of the v2 pipeline, collected by
collect_identity_full_site_tables.py), and the reference .fai for chromosome lengths."""
import collections
from fig_common import *   # Phase N: print size, >= 7 pt text, RGB TIFF via fig_common.save (also honours FIG_DIR)
BIN=20_000_000; CHR=[str(i) for i in range(1,20)]
LEN={l.split("\t")[0]:int(l.split("\t")[1]) for l in open(f"{PROJ}/ref/b6_fasta/GRCm39.primary.fa.fai") if l.split("\t")[0] in CHR}
T=f"{PROJ}/analysis_v2/results/identity/full_site_tables"; cand=pd.read_csv(f"{T}/candidates.tsv", sep="\t"); ref=cand[cand.chosen=="Y"].set_index("strain")
TL={"adipose":"mesenteric adipose","hypothalamus":"hypothalamus","liver":"liver"}
rows=[(f"GSE55788 'C57BL/6J' heart {i} (GSM{1346086+i})",f"{T}/heart_B6_C57_M_12wk_mRNA_{i}.ase.tsv") for i in range(1,5)]+\
     [(f"Reference {st} ({TL[ref.tissue[st]]}, {ref.GSE[st]}, {ref.GSM[st]})",f"{T}/{ref.sample_id[st]}.ase.tsv") for st in ("C57BL/6J","FVB/NJ")]
cols=[(c,b) for c in CHR for b in range(int(np.ceil(LEN[c]/BIN)))]
mat=np.full((len(rows),len(cols)),np.nan); idx={k:i for i,k in enumerate(cols)}
for i,(lab,p) in enumerate(rows):
    a=pd.read_csv(p, sep="\t", dtype={"contig":str}); a=a[(a.refCount+a.altCount)>=20]; a=a[a.contig.isin(CHR)]
    a["af"]=a.altCount/(a.refCount+a.altCount); a["bin"]=a.position//BIN
    g=a.groupby(["contig","bin"]).af.agg(n="size", hi=lambda x:(x>0.9).mean())
    for (c,b),r in g.iterrows():
        if r.n>=10: mat[i,idx[(c,int(b))]]=r.hi
fig,ax=plt.subplots(figsize=(FULL_W,3.3)); fig.subplots_adjust(left=0.3,right=0.88,top=0.8,bottom=0.15)
im=ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0].replace(" (","\n(",1).replace(", GSE",",\nGSE",1) for r in rows], fontsize=7, linespacing=0.95)   # display only; the _data.tsv index keeps the one-line labels
starts=[]; pos=0
for c in CHR:
    n=int(np.ceil(LEN[c]/BIN)); starts.append((pos+n/2-0.5,c)); pos+=n
    if c!="19": ax.axvline(pos-0.5, color="white", lw=1.2)
ax.set_xticks([s for s,_ in starts]); ax.set_xticklabels([c for _,c in starts], fontsize=7); ax.set_xlabel("Chromosome (20-Mb bins)", fontsize=7.5)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
cb=fig.colorbar(im, ax=ax, fraction=0.03, pad=0.015); cb.set_label("Fraction of FVB/NJ-informative\nsites homozygous for the\nFVB/NJ allele", fontsize=7); cb.ax.tick_params(labelsize=7)
ax.set_title("Homozygous FVB/NJ-allele blocks in GSE55788 \"wild-type C57BL/6J\" heart RNA-seq\n(bins with ≥ 10 sites at depth ≥ 20)", fontsize=8)
save(fig, "SFig1_GSE55792_haplotype_blocks", pd.DataFrame(mat, index=[r[0] for r in rows], columns=[f"chr{c}:{b*20}-{(b+1)*20}Mb" for c,b in cols]))
