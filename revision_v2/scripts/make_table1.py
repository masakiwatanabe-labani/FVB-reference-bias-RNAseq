#!/usr/bin/env python3
"""Table 1 generated from revision_v2/manifest/manifest.tsv (+ analysis_v2/results/qc/sample_identity_verification.tsv when present).
Output: manuscript/revised_v2/generated/table1.md  (rule 0.4-1: never hand-edited; Phase F re-generates and diffs)."""
import csv, collections, os
PROJ="/usr/local/jupyter/FVB_bias"; OUT=f"{PROJ}/manuscript/revised_v2/generated/table1.md"
M=[r for r in csv.DictReader(open(f"{PROJ}/revision_v2/manifest/manifest.tsv"), delimiter="\t") if r["include"]=="Y"]
qc={}
p=f"{PROJ}/analysis_v2/results/qc/sample_identity_verification.tsv"
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter="\t"): qc[r["sample_id"]]=r
TISSUE_LABEL={"adipose":"Mesenteric adipose","hypothalamus":"Hypothalamus","liver":"Liver","skin":"Back skin","mammary_tumour_PyVT":"Mammary tumour (MMTV-PyVT)","keratinocytes":"Primary keratinocytes (LSL-KrasG12D, Cre-transduced)","cerebral_cortex":"Cerebral (frontal) cortex"}
COHORT={"GSE123875":"Fructose study (PRJNA510125; SubSeries of GSE123896)","GSE123893":"Fructose study (PRJNA510129; SubSeries of GSE123896)","GSE123894":"Fructose study (PRJNA510133; SubSeries of GSE123896)","GSE86315":"Imiquimod strain panel (PRJNA341396)"}
def nq(rows):
    if not qc: return "pending"
    ok=sum(1 for r in rows if qc.get(r["sample_id"],{}).get("match")=="Y"); return f"{ok}/{len(rows)}"
def cond(rows):
    c=collections.Counter((r["diet"], r["sex"]) for r in rows)
    parts=[]
    for (d,s),n in sorted(c.items()): parts.append(f"{d}{(' '+s) if s else ''} ×{n}")
    return "; ".join(parts)
L=["# Table 1. RNA-seq datasets used in this study","","Generated from the sample manifest (`revision_v2/manifest/manifest.tsv`) by `revision_v2/scripts/make_table1.py`; *n* = biological animals (technical sequencing runs of one animal were merged).","",
   "## Three-arm strain comparison datasets (FVB/NJ versus C57BL/6J)","",
   "| GEO accession | Tissue | Cohort | *n* FVB/NJ (conditions) | *n* C57BL/6J (conditions) | Library | Read length | Runs per animal | Sequence identity QC passed |","|---|---|---|---|---|---|---|---|---|"]
sc=[r for r in M if r["role"]=="strain_comparison"]
for gse in sorted(set(r["GSE"] for r in sc), key=lambda g:[r["tissue"] for r in sc if r["GSE"]==g][0]):
    rows=[r for r in sc if r["GSE"]==gse]; fv=[r for r in rows if r["strain"]=="FVB/NJ"]; b6=[r for r in rows if r["strain"]=="C57BL/6J"]
    runs=sorted(set(r["n_technical_runs"] for r in rows)); rl=sorted(set(str(int(float(r['read_length']))) for r in rows if r['read_length'])); lay=sorted(set(r["library_layout"] for r in rows))
    L.append(f"| {gse} | {TISSUE_LABEL[rows[0]['tissue']]} | {COHORT.get(gse,'')} | {len(fv)} ({cond(fv)}) | {len(b6)} ({cond(b6)}) | {'/'.join(x.capitalize()+'-end' for x in lay)} | {'/'.join(rl)} bp | {'/'.join(runs)} | {nq(rows)} |")
L+=["","## Negative-control dataset (no FVB/NJ-versus-C57BL/6J contrast)","","| GEO accession | Tissue | Strain (GEO attribute) | *n* | Conditions | Library | Read length | Sequence identity QC |","|---|---|---|---|---|---|---|---|"]
nc=[r for r in M if r["role"]=="negative_control_candidate"]
if nc:
    L.append(f"| GSE175625 | {TISSUE_LABEL['mammary_tumour_PyVT']} | B6.FVB-Tg(MMTV-PyVT)634Mul/LellJ (C57BL/6 congenic) | {len(nc)} | {cond(nc)} | Single-end | {int(float(nc[0]['read_length']))} bp | {nq(nc)} (expected: C57BL/6) |")
L+=["","## F1 hybrid datasets used for allele-specific expression analysis (Fig. 2)","","| GEO accession | Tissue | Cross (GEO attribute) | *n* (F1) | Library | Read length | Sequence identity QC |","|---|---|---|---|---|---|---|"]
for gse in ("GSE293440","GSE200632"):
    rows=[r for r in M if r["role"]=="f1_ase" and r["GSE"]==gse]
    if rows: L.append(f"| {gse} | {TISSUE_LABEL[rows[0]['tissue']]} | {rows[0]['strain']} | {len(rows)} | {rows[0]['library_layout'].capitalize()}-end | {int(float(rows[0]['read_length']))} bp | {nq(rows)} (expected: F1) |")
L+=["","Sequence identity QC: allele fractions at FVB/NJ-homozygous exonic SNP sites (203,510; a fixed 10% subset for all but the F1 animals) in the B6-reference alignment (see *Sample identity verification*); full per-sample results in Supplementary Table S1.",
    "Sex and age are not recorded in the GEO/BioSample attributes of the fructose study; the imiquimod study records sex (M/F). All FVB/NJ animals are from the Jackson Laboratory substrain, matched to the FVB_NJ variant set."]
open(OUT,"w").write("\n".join(L)+"\n"); print(open(OUT).read())
