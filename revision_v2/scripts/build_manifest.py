#!/usr/bin/env python3
"""Phase B-2: build revision_v2/manifest/manifest.tsv — one row per biological sample (BioSample).
Inputs: revision_v2/manifest/fructose_all_runs.tsv (B-1), revision_v2/audit/old_runs_metadata.tsv (A-1),
        analysis_v2/qc/strain_identity_qc.tsv (A-2, old BAMs only).
Selection rules are code, not hand edits:
  strain comparison: GEO sample attribute strain == "FVB/NJ" or "C57BL/6J" (exact), tissue from the series,
                     design_main = all animals of the tissue (primary model, both diets/treatments); design_sensitivity =
                     control-diet/control-treatment animals only (control-only sensitivity model). Phase U-2 (2026-10-24): the two
                     flags were inverted relative to the Methods and to the models fitted (bias_analysis.py: "primary" on all 70
                     animals, "water" on the 35 control animals); swapped here, and the include_reason texts of the
                     water and fructose rows reworded to match (IMQ and control-cream texts were already correct).
  DBA/2J: excluded. GSE135230: excluded (single mixed background). GSE175625: candidate negative control
  (include = PENDING until the human decides, rule B-5). F1 (GSE200632, GSE293440): include = Y (reused).
"""
import collections, csv
PROJ="/usr/local/jupyter/FVB_bias"; MAN=f"{PROJ}/revision_v2/manifest"; AUD=f"{PROJ}/revision_v2/audit"
TISSUE={"GSE123875":"adipose","GSE123893":"hypothalamus","GSE123894":"liver"}
STRAIN_CODE={"FVB/NJ":"FVB","C57BL/6J":"B6","DBA/2J":"DBA"}
qc={}
for r in csv.DictReader(open(f"{PROJ}/analysis_v2/qc/strain_identity_qc.tsv"), delimiter="\t"): qc[r["SRR"]]=r["qc_verdict"]

def agg_qc(srrs):
    v=sorted({qc[s] for s in srrs if s in qc})
    return ";".join(v) if v else ""

samples={}  # BioSample -> dict
# --- fructose study (all runs)
for r in csv.DictReader(open(f"{MAN}/fructose_all_runs.tsv"), delimiter="\t"):
    s=samples.setdefault(r["BioSample"], dict(GSE=r["GSE"], GSM=r["GSM"], BioSample=r["BioSample"], SRR=[], strain=r["geo_strain"],
        tissue=TISSUE[r["GSE"]], diet=r["geo_treatment"].lower(), sex=r["biosample_sex"], age=r["biosample_age"],
        library_layout=r["library_layout"], read_length=r["read_length_est"], title=r["sample_title"], bytes=0, reads=0, ftp=[], md5=[], role="strain_comparison"))
    s["SRR"].append(r["SRR"]); s["bytes"]+=sum(int(x) for x in r["fastq_bytes"].split(";") if x); s["reads"]+=int(r["read_count"] or 0)
    s["ftp"].append(r["fastq_ftp"]); s["md5"].append(r["fastq_md5"])
# --- old-analysis samples not in the fructose study (F1, PyVT, microglia); fastq fields from the ENA study run lists saved in Phase A
import glob
enarun={}
for srp in ("SRP217051","SRP321579","SRP369603","SRP573960"):
    for x in csv.DictReader(open(sorted(glob.glob(f"{AUD}/raw/ena_study_runs_{srp}_*.tsv"))[-1]), delimiter="\t"): enarun[x["run_accession"]]=x
for r in csv.DictReader(open(f"{AUD}/old_runs_metadata.tsv"), delimiter="\t"):
    if r["GSE"] in TISSUE: continue
    tissue={"GSE135230":"microglia","GSE175625":"mammary_tumour_PyVT","GSE200632":"keratinocytes","GSE293440":"cerebral_cortex"}[r["GSE"]]
    strain=r["geo_strain"] or r["biosample_strain"]
    s=samples.setdefault(r["BioSample"], dict(GSE=r["GSE"], GSM=r["GSM"], BioSample=r["BioSample"], SRR=[], strain=strain, tissue=tissue,
        diet=r["geo_treatment"] or r["geo_genotype"], sex=r["geo_sex"], age=r["geo_age"], library_layout=r["library_layout"], read_length="",
        title=r["geo_title"], bytes=0, reads=0, ftp=[], md5=[], role={"GSE135230":"excluded","GSE175625":"negative_control_candidate"}.get(r["GSE"],"f1_ase")))
    s["SRR"].append(r["SRR"]); s["reads"]+=int(r["read_count"] or 0)
    e=enarun.get(r["SRR"])
    if e:
        s["bytes"]+=sum(int(x) for x in e["fastq_bytes"].split(";") if x); s["ftp"].append(e["fastq_ftp"]); s["md5"].append(e["fastq_md5"])
        rc=int(e["read_count"] or 0); bc=int(e["base_count"] or 0); s["read_length"]=round(bc/rc/(2 if e["library_layout"]=="PAIRED" else 1),1) if rc else ""

# --- GSE55792 heart (added after Phase B STOP; decision 3 in design_proposal.md): bulk polyA whole-heart RNA-seq,
#     GEO strain attribute C57BL/6J (GSE55788, 12 wk males) or FVB/NJ (GSE55790, males 6-8 wk = closest ages to 12 wk).
#     Source: revision_v2/audit/gse55792_runs.tsv (GEO SOFT + ENA, Phase B). Single-end 50 bp, HiSeq 2000, 1 run per animal.
for r in csv.DictReader(open(f"{AUD}/gse55792_runs.tsv"), delimiter="\t"):
    ch={x.partition(": ")[0].strip().lower():x.partition(": ")[2].strip() for x in r["characteristics"].split("|")}
    if r["library_strategy"]!="RNA-Seq" or r["library_selection"]!="cDNA" or "Whole-heart" not in ch.get("tissue",""): continue
    strain=ch.get("strain",""); sex=ch.get("sex") or ch.get("gender",""); age=ch.get("age (weeks)","")
    s=samples.setdefault(r["sample_accession"], dict(GSE="GSE55792", GSM=r["GSM"], BioSample=r["sample_accession"], SRR=[], strain=strain, tissue="heart",
        diet="standard chow", sex=sex, age=f"{age} wk", library_layout=r["library_layout"], read_length="", title=r["title"], bytes=0, reads=0, ftp=[], md5=[], role="strain_comparison"))
    s["SRR"].append(r["run_accession"]); s["reads"]+=int(r["read_count"] or 0); s["bytes"]+=sum(int(x) for x in r["fastq_bytes"].split(";") if x)
    s["ftp"].append(r["fastq_ftp"]); s["md5"].append(r["fastq_md5"])
    rc=int(r["read_count"] or 0); bc=int(r["base_count"] or 0); s["read_length"]=round(bc/rc) if rc else ""

# --- GSE86315 skin (reserve candidate adopted 2026-10-14 after GSE55792 failed QC; candidate_datasets.md): back skin, 7 strains x sex x
#     control/IMQ, 2 animals per cell. Selected: GEO strain attribute FVB/NJ or C57BL/6J (16 animals: 2 sexes x 2 treatments x 2).
#     Source: revision_v2/audit/gse86315_runs.tsv. Single-end 52 bp GA IIx, 1 run per animal. "diet" column holds the treatment.
for r in csv.DictReader(open(f"{AUD}/gse86315_runs.tsv"), delimiter="\t"):
    ch={x.partition(": ")[0].strip().lower():x.partition(": ")[2].strip() for x in r["characteristics"].split("|")}
    strain=ch.get("strain","")
    if strain not in ("FVB/NJ","C57BL/6J"): continue
    s=samples.setdefault(r["sample_accession"], dict(GSE="GSE86315", GSM=r["GSM"], BioSample=r["sample_accession"], SRR=[], strain=strain, tissue="skin",
        diet=ch.get("treatment",""), sex={"male":"M","female":"F"}.get(ch.get("sex","").lower(),ch.get("sex","")), age="", library_layout=r["library_layout"], read_length="", title=r["title"], bytes=0, reads=0, ftp=[], md5=[], role="strain_comparison"))
    s["SRR"].append(r["run_accession"]); s["reads"]+=int(r["read_count"] or 0); s["bytes"]+=sum(int(x) for x in r["fastq_bytes"].split(";") if x)
    s["ftp"].append(r["fastq_ftp"]); s["md5"].append(r["fastq_md5"])
    rc=int(r["read_count"] or 0); bc=int(r["base_count"] or 0); s["read_length"]=round(bc/rc) if rc else ""

def decide(s):
    """returns include, include_reason, design_main, design_sensitivity"""
    if s["role"]=="strain_comparison" and s["tissue"]=="heart":
        # 2026-10-14: all four GSE55788 "C57BL/6J" animals failed the sequence identity QC (15-18 % homozygous FVB-allele sites in
        # large chromosomal blocks = mixed background; revision_v2/audit/gse55792_qc_note.md) -> no B6 comparator -> dataset not used.
        return "N","GSE55792 not used: the C57BL/6J group (GSE55788, 4/4 animals) failed the sequence-based strain identity QC (mixed background); see audit/gse55792_qc_note.md","N","N"
        if s["strain"] not in ("FVB/NJ","C57BL/6J"): return "N",f"strain '{s['strain']}'","N","N"
        if s["sex"]!="M": return "N","heart: female animals not used (B6 reference group is male only)","N","N"
        if s["strain"]=="FVB/NJ" and s["age"] not in ("6 wk","7 wk","8 wk"): return "N",f"heart FVB/NJ: age {s['age']} outside 6-8 wk window (closest to the 12 wk B6 group)","N","N"
        return "Y",f"{s['strain']} whole heart male {s['age']}: additional tissue (no diet factor; model ~ strain)","Y","Y"
    if s["role"]=="strain_comparison" and s["tissue"]=="skin":
        if s["diet"]=="control": return "Y",f"{s['strain']} skin {s['sex']} control cream: additional tissue, main design (~ treatment + sex + strain) and control-only sensitivity","Y","Y"
        return "Y",f"{s['strain']} skin {s['sex']} IMQ-treated: additional tissue, main design only","Y","N"
    if s["role"]=="strain_comparison":
        if s["strain"]=="DBA/2J": return "N","strain DBA/2J: not part of the FVB/NJ vs C57BL/6J comparison","N","N"
        if s["strain"] not in ("FVB/NJ","C57BL/6J"): return "N",f"strain '{s['strain']}' not FVB/NJ or C57BL/6J","N","N"
        if s["diet"]=="water": return "Y",f"{s['strain']} {s['tissue']} control diet (water): main design and control-only sensitivity","Y","Y"
        if s["diet"]=="fructose": return "Y",f"{s['strain']} {s['tissue']} fructose diet: main design (~ diet + strain) only","Y","N"
        return "N",f"unexpected diet '{s['diet']}'","N","N"
    if s["role"]=="f1_ase": return "Y","F1 hybrid for allele-specific expression (Phase A QC: F1)","-","-"
    if s["role"]=="negative_control_candidate": return "Y","negative control (decision 2026-10-14, Option A): B6.FVB-Tg(MMTV-PyVT) congenic tumours, Phase A QC 'B6 pure' ×13; evaluated with the per-sample reference-swap ratio, not the strain-contrast bias","-","-"
    return "N","GSE135230: all animals on one mixed B6.129 × FVB/N-Tg background (QC F1-like / backcross-like); no FVB vs B6 contrast","-","-"

cols=["sample_id","GSE","GSM","BioSample","SRR","strain","tissue","diet","sex","age","library_layout","read_length","include","include_reason","qc_verdict",
      "role","design_main","design_sensitivity","n_technical_runs","read_count_total","fastq_bytes_total","fastq_ftp","fastq_md5"]
with open(f"{MAN}/manifest.tsv","w") as o:
    o.write("\t".join(cols)+"\n")
    for s in sorted(samples.values(), key=lambda s:(s["role"]!="strain_comparison", s["GSE"], s["tissue"], s["strain"], s["diet"], s["title"])):
        inc,reason,dm,ds=decide(s)
        code=STRAIN_CODE.get(s["strain"],"X")
        sid=f"{s['tissue']}_{code}_{s['title'].replace(' ','').replace(',','_')}"
        o.write("\t".join([sid,s["GSE"],s["GSM"],s["BioSample"],",".join(sorted(s["SRR"])),s["strain"],s["tissue"],s["diet"],s["sex"],s["age"],
            s["library_layout"],str(s["read_length"]),inc,reason,agg_qc(s["SRR"]),s["role"],dm,ds,str(len(s["SRR"])),str(s["reads"]),str(s["bytes"]),
            ";".join(x for x in s["ftp"] if x),";".join(x for x in s["md5"] if x)])+"\n")
# summary
c=collections.Counter()
for s in samples.values():
    inc,_,dm,ds=decide(s)
    if s["role"]=="strain_comparison" and inc=="Y": c[(s["tissue"],s["strain"],s["diet"])]+=1
for k,v in sorted(c.items()): print(k,v)
print("manifest rows:",len(samples))
