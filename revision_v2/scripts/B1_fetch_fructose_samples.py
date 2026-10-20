#!/usr/bin/env python3
"""Phase B-1: enumerate every run of the three fructose-study GEO series (GSE123875 adipose, GSE123893
hypothalamus, GSE123894 liver) from the ENA study run lists saved in Phase A, fetch the BioSample XML
attributes for every sample (dated raw copies in revision_v2/audit/raw/), join with GEO SOFT
characteristics, and write revision_v2/manifest/fructose_all_runs.tsv (one row per run).
"""
import os, glob, gzip, time, re, datetime, collections, urllib.request
from xml.etree import ElementTree as ET
PROJ="/usr/local/jupyter/FVB_bias"; AUD=f"{PROJ}/revision_v2/audit"; RAW=f"{AUD}/raw"; MAN=f"{PROJ}/revision_v2/manifest"
DATE=datetime.date.today().isoformat(); os.makedirs(MAN, exist_ok=True)
STUDIES={"SRP173496":"GSE123875","SRP173500":"GSE123893","SRP173501":"GSE123894"}

def fetch(url, retries=4):
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r: return r.read()
        except Exception as e:
            print("  retry", i+1, url, e, flush=True); time.sleep(2*(i+1))
    raise RuntimeError(url)

# GEO SOFT samples (all series samples)
geo={}; by_biosample={}
for gse in STUDIES.values():
    p=sorted(glob.glob(f"{RAW}/{gse}_family.soft_*.gz"))[-1]; cur=None
    with gzip.open(p,"rt",errors="replace") as f:
        for line in f:
            line=line.rstrip("\n")
            if line.startswith("^SAMPLE"): cur=line.split("=")[1].strip(); geo[cur]=collections.defaultdict(list); geo[cur]["GSE"]=[gse]
            elif cur and line.startswith("!Sample_"):
                k,_,v=line[1:].partition(" = "); geo[cur][k].append(v)
for gsm,d in geo.items():
    for rel in d.get("Sample_relation",[]):
        m=re.search(r"(SAMN\d+)",rel)
        if m: by_biosample[m.group(1)]=gsm
def ch(d,key):
    for c in d.get("Sample_characteristics_ch1",[]):
        k,_,v=c.partition(": ")
        if k.strip().lower()==key: return v.strip()
    return ""

rows=[]
for srp,gse in STUDIES.items():
    p=sorted(glob.glob(f"{RAW}/ena_study_runs_{srp}_*.tsv"))[-1]
    lines=open(p).read().splitlines(); h=lines[0].split("\t")
    for l in lines[1:]:
        r=dict(zip(h,l.split("\t"))); r["GSE"]=gse; rows.append(r)
print("runs:",len(rows))

# BioSample attributes (reuse Phase A XMLs when present)
attrs={}
for s in sorted({r["sample_accession"] for r in rows}):
    xmls=glob.glob(f"{RAW}/biosample_{s}_*.xml")
    if xmls: xml=open(xmls[-1],"rb").read()
    else:
        xml=fetch(f"https://www.ebi.ac.uk/ena/browser/api/xml/{s}"); open(f"{RAW}/biosample_{s}_{DATE}.xml","wb").write(xml); time.sleep(1)
    root=ET.fromstring(xml); a={}
    for smp in root.iter("SAMPLE"):
        a["_alias"]=smp.attrib.get("alias",""); a["_title"]=(smp.findtext("TITLE") or "").strip()
        for x in smp.iter("SAMPLE_ATTRIBUTE"):
            a[(x.findtext("TAG") or "").strip().lower()]=(x.findtext("VALUE") or "").strip()
    attrs[s]=a
print("samples:",len(attrs))

cols=["GSE","GSM","BioSample","SRS","SRR","SRX","sample_title","geo_strain","geo_treatment","geo_source","biosample_strain","biosample_treatment",
      "biosample_sex","biosample_age","library_layout","instrument","read_count","base_count","read_length_est","fastq_bytes","fastq_ftp","fastq_md5"]
with open(f"{MAN}/fructose_all_runs.tsv","w") as o:
    o.write("\t".join(cols)+"\n")
    for r in sorted(rows,key=lambda r:(r["GSE"],r["sample_accession"],r["run_accession"])):
        s=r["sample_accession"]; gsm=by_biosample.get(s,""); g=geo.get(gsm,{}); a=attrs[s]
        rc=int(r["read_count"] or 0); bc=int(r["base_count"] or 0); npair=2 if r["library_layout"]=="PAIRED" else 1
        rl=round(bc/rc/npair,1) if rc else ""
        o.write("\t".join([r["GSE"],gsm,s,r["secondary_sample_accession"],r["run_accession"],r["experiment_accession"],
            " ".join(g.get("Sample_title",[])) or r["sample_title"], ch(g,"strain"), ch(g,"treatment"), " ".join(g.get("Sample_source_name_ch1",[])),
            a.get("strain",""), a.get("treatment",""), a.get("sex",""), a.get("age",""),
            r["library_layout"], r["instrument_model"], r["read_count"], r["base_count"], str(rl), r["fastq_bytes"], r["fastq_ftp"], r["fastq_md5"]])+"\n")
print("B1 DONE")
