#!/usr/bin/env python3
"""Merge per-sample featureCounts outputs into per-tissue matrices and build the sample QC / mapping tables.
Inputs: revision_v2/manifest/manifest.tsv (include=Y rows that have done.all) and analysis_v2/data/<tissue>/<sample>/.
Outputs (analysis_v2/results/):
  counts/<tissue>.<arm>.<mode>.tsv        gene x sample raw counts (arm: naive|wasp|fvb; mode: unique|multi)
  counts/<tissue>.samples.tsv             sample sheet (sample_id, strain, diet, sex, age, layout, read_length)
  qc/sample_identity_verification.tsv     Supplementary Table "Sample identity verification": all include=Y samples (metrics on the fixed 10 % site
                                          subset, qc/identity.sub10pct.tsv where the sample was run before speed-up Step 3, else qc/identity.tsv),
                                          plus every processed sample that FAILED the QC (include=N with a qc/identity*.tsv; GSE55792 heart), plus
                                          the F1 samples re-checked in Phase H-1 (full site set). Column site_set records which site set the metrics come from.
  qc/mapping_stats_per_sample.tsv         STAR / WASP / featureCounts summary per sample and arm
"""
import csv, os, re, json, collections
PROJ="/usr/local/jupyter/FVB_bias"; DATA=f"{PROJ}/analysis_v2/data"
RES=os.environ.get("RES_DIR", f"{PROJ}/analysis_v2/results")   # RES_DIR: per-tissue dry runs write to a scratch dir; analysis_v2/results/ is written only by the final run
os.makedirs(f"{RES}/counts", exist_ok=True); os.makedirs(f"{RES}/qc", exist_ok=True)
M=[r for r in csv.DictReader(open(f"{PROJ}/revision_v2/manifest/manifest.tsv"), delimiter="\t")]
def read_counts(p):
    d={}
    with open(p) as f:
        for l in f:
            if l.startswith("#") or l.startswith("Geneid"): continue
            r=l.rstrip("\n").split("\t"); d[r[0]]=r[6]
    return d
def star_stats(p):
    s={}
    for l in open(p):
        if "|" in l:
            k,v=[x.strip() for x in l.split("|",1)]; s[k]=v
    return s
def fc_summary(p):
    s={}
    for l in open(p):
        r=l.rstrip("\n").split("\t")
        if len(r)>=2 and r[0]!="Status": s[r[0]]=int(r[1])
    return s
def identity(d):
    """(last row of the identity table, site_set). Prefer the 10 % subset table so that every sample is evaluated on the same
    20,351 sites (samples run before 2026-10-16 have the full-set identity.tsv plus identity.sub10pct.tsv; later ones only identity.tsv on the subset)."""
    p=f"{d}/qc/identity.sub10pct.tsv" if os.path.exists(f"{d}/qc/identity.sub10pct.tsv") else f"{d}/qc/identity.tsv"
    row=open(p).read().splitlines()[-1].split("\t")
    # 2026-10-20: the label must follow the file actually used. hypothalamus_B6_CWHyp3 (10th animal, run with the full set) had no
    # subset table and its 33,355-site row was labelled as the subset; its identity.sub10pct.tsv was re-derived from the full ase.tsv
    # (same verdict) and a full-set table (> 20,351 possible sites) now stops the merge instead of being mislabelled.
    if not p.endswith("identity.sub10pct.tsv") and int(row[0])>20351: raise SystemExit(f"identity(): {p} has {row[0]} sites (full set) and no qc/identity.sub10pct.tsv; re-derive the subset table first")
    return row, "sub10pct (20,351 sites)"
done=[r for r in M if r["include"]=="Y" and os.path.exists(f"{DATA}/{r['tissue']}/{r['sample_id']}/done.all")]
by_t=collections.defaultdict(list)
for r in done: by_t[r["tissue"]].append(r)
qc_rows=[]; map_rows=[]
for t,rows in by_t.items():
    with open(f"{RES}/counts/{t}.samples.tsv","w") as o:
        o.write("sample_id\tstrain\tdiet\tsex\tage\tlibrary_layout\tread_length\tGSE\tGSM\n")
        for r in rows: o.write("\t".join([r["sample_id"],r["strain"],r["diet"],r["sex"],r["age"],r["library_layout"],r["read_length"],r["GSE"],r["GSM"]])+"\n")
    for arm,mode in [("naive","unique"),("wasp","unique"),("fvb","unique"),("naive","multi"),("fvb","multi")]:
        mats=[]; ids=[]
        for r in rows:
            p=f"{DATA}/{t}/{r['sample_id']}/counts/{r['sample_id']}.{arm}.{mode}.txt"
            mats.append(read_counts(p)); ids.append(r["sample_id"])
        genes=list(mats[0].keys())
        with open(f"{RES}/counts/{t}.{arm}.{mode}.tsv","w") as o:
            o.write("gene_id\t"+"\t".join(ids)+"\n")
            for g in genes: o.write(g+"\t"+"\t".join(m[g] for m in mats)+"\n")
    for r in rows:
        d=f"{DATA}/{t}/{r['sample_id']}"; sid=r["sample_id"]
        q,ss=identity(d)
        fp=json.load(open(f"{d}/qc/fastp.json"))["summary"]
        qc_rows.append([sid,r["GSE"],r["GSM"],r["BioSample"],r["SRR"],t,r["strain"],r["diet"],r["sex"],r["age"],r["library_layout"],
                        str(fp["before_filtering"]["total_reads"]),str(fp["after_filtering"]["total_reads"])]+q[:6]+[q[6],q[7],ss,"Y",""])
        for arm in ("naive","fvb"):
            s=star_stats(f"{d}/bam/{sid}.{arm}.Log.final.out")
            fcu=fc_summary(f"{d}/counts/{sid}.{arm}.unique.txt.summary"); fcm=fc_summary(f"{d}/counts/{sid}.{arm}.multi.txt.summary")
            map_rows.append([sid,t,r["strain"],r["diet"],arm,s["Number of input reads"],s["Uniquely mapped reads number"],s["Uniquely mapped reads %"],
                             s["% of reads mapped to multiple loci"],s["% of reads unmapped: too short"],str(fcu.get("Assigned","")),str(fcu.get("Unassigned_MultiMapping","")),str(round(fcm.get("Assigned",0)))])
        w=fc_summary(f"{d}/counts/{sid}.wasp.unique.txt.summary")
        vw={l.split("\t")[0]:l.split("\t")[1] for l in open(f"{d}/qc/vW_tag_counts.tsv").read().splitlines()}
        map_rows.append([sid,t,r["strain"],r["diet"],"wasp","","","","","",str(w.get("Assigned","")),str(w.get("Unassigned_MultiMapping","")),"vW1="+vw.get("1","0")+";vW2="+vw.get("2","0")+";vW3-7="+str(sum(int(v) for k,v in vw.items() if k not in ("1","2","none")))+";none="+vw.get("none","0")])
# processed samples that FAILED the identity QC (rule 0.4-4: failing samples are not used but are reported): include=N rows with a QC table
for r in M:
    d=f"{DATA}/{r['tissue']}/{r['sample_id']}"
    if r["include"]=="N" and os.path.exists(f"{d}/qc/identity.tsv"):
        q,ss=identity(d); fp=json.load(open(f"{d}/qc/fastp.json"))["summary"]
        qc_rows.append([r["sample_id"],r["GSE"],r["GSM"],r["BioSample"],r["SRR"],r["tissue"],r["strain"],r["diet"],r["sex"],r["age"],r["library_layout"],
                        str(fp["before_filtering"]["total_reads"]),str(fp["after_filtering"]["total_reads"])]+q[:6]+[q[6],q[7],ss,"N",r["include_reason"]])
# F1 rows (Phase H-1, 2026-10-21: identity re-check on the re-aligned F1 BAMs, full 203,510-site set, written by f1_ase_rerun.sh to
# analysis_v2/results/f1_ase/v2/<sid>.identity.tsv; read counts from the new fastp run). Replaces the Phase A rows (old BAMs deleted;
# the Phase A rows are kept in revision_v2/audit/f1_old/sample_identity_verification.phaseA_F1rows.tsv).
for r in M:
    if r["role"]=="f1_ase" and r["include"]=="Y":
        sid=r["sample_id"]; d=f"{DATA}/{r['tissue']}/{sid}"
        q=open(f"{PROJ}/analysis_v2/results/f1_ase/v2/{sid}.identity.tsv").read().splitlines()[-1].split("\t")
        fp=json.load(open(f"{d}/qc/fastp.json"))["summary"]
        qc_rows.append([sid,r["GSE"],r["GSM"],r["BioSample"],r["SRR"],r["tissue"],r["strain"],r["diet"],r["sex"],r["age"],r["library_layout"],
                        str(fp["before_filtering"]["total_reads"]),str(fp["after_filtering"]["total_reads"]),
                        q[0],q[1],q[2],q[3],q[4],q[5],"F1 (B6 x FVB)","Y" if q[5]=="F1 (B6 x FVB)" else "N","full (203,510 sites)","Y",""])
with open(f"{RES}/qc/sample_identity_verification.tsv","w") as o:
    o.write("\t".join(["sample_id","GSE","GSM","BioSample","SRR","tissue","strain_GEO","diet_or_group","sex","age","library_layout","reads_raw","reads_after_fastp",
                       "n_sites","mean_alt_fraction","frac_hi","frac_mid","frac_lo","qc_verdict","expected","match","site_set","used_in_analysis","exclusion_reason"])+"\n")
    for r in qc_rows: o.write("\t".join(r)+"\n")
with open(f"{RES}/qc/mapping_stats_per_sample.tsv","w") as o:
    o.write("\t".join(["sample_id","tissue","strain","diet","arm","input_reads","uniq_mapped","uniq_pct","multi_pct","unmapped_too_short_pct","fc_assigned_unique","fc_unassigned_multimapping","fc_assigned_multi_fraction_or_vW"])+"\n")
    for r in map_rows: o.write("\t".join(r)+"\n")
print("tissues:", {t:len(v) for t,v in by_t.items()}, "| qc rows:", len(qc_rows))
