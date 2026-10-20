#!/usr/bin/env python3
"""Phase A-1: fetch ENA run reports, ENA/BioSample XML attributes and GEO SOFT files
for every SRR used in the old analysis. All raw responses are saved (dated) under
revision_v2/audit/raw/. curl is not installed on this machine, so urllib is used.
"""
import sys, time, gzip, re, json, os, datetime, urllib.request, urllib.error
from xml.etree import ElementTree as ET

PROJ = "/usr/local/jupyter/FVB_bias"
AUD = f"{PROJ}/revision_v2/audit"
RAW = f"{AUD}/raw"
DATE = datetime.date.today().isoformat()
os.makedirs(RAW, exist_ok=True)

def fetch(url, retries=4, sleep=2):
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception as e:
            print(f"  retry {i+1} {url}: {e}", flush=True)
            time.sleep(sleep * (i + 1))
    raise RuntimeError(f"failed: {url}")

srrs = [l.strip() for l in open(f"{AUD}/old_srr_list.txt") if l.strip()]
FIELDS = ("run_accession,sample_accession,secondary_sample_accession,sample_alias,sample_title,"
          "experiment_accession,study_accession,secondary_study_accession,library_layout,"
          "library_strategy,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes")

# 1) per-SRR ENA filereport
out = open(f"{AUD}/old_runs_ena.tsv", "w")
header_written = False
for srr in srrs:
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={srr}&result=read_run&fields={FIELDS}&format=tsv"
    txt = fetch(url).decode()
    open(f"{RAW}/ena_filereport_{srr}_{DATE}.tsv", "w").write(txt)
    lines = [l for l in txt.splitlines() if l.strip()]
    if not header_written:
        out.write(lines[0] + "\n"); header_written = True
    for l in lines[1:]:
        out.write(l + "\n")
    print(f"[ENA run] {srr} ok", flush=True)
    time.sleep(1)
out.close()

# 2) per-sample ENA XML (sample_accession = SAMN / secondary = SRS)
rows = [l.rstrip("\n").split("\t") for l in open(f"{AUD}/old_runs_ena.tsv")]
hdr = rows[0]; rows = rows[1:]
idx = {h: i for i, h in enumerate(hdr)}
samples = sorted({r[idx["sample_accession"]] for r in rows})
attr_out = open(f"{AUD}/old_samples_attributes.tsv", "w")
attr_out.write("sample_accession\tsecondary_sample_accession\talias\ttitle\ttag\tvalue\n")
for s in samples:
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{s}"
    xml = fetch(url)
    open(f"{RAW}/biosample_{s}_{DATE}.xml", "wb").write(xml)
    root = ET.fromstring(xml)
    for smp in root.iter("SAMPLE"):
        alias = smp.attrib.get("alias", "")
        title = (smp.findtext("TITLE") or "").strip()
        sec = ""
        for ident in smp.iter("SECONDARY_ID"):
            sec = ident.text or ""
        for a in smp.iter("SAMPLE_ATTRIBUTE"):
            tag = (a.findtext("TAG") or "").strip(); val = (a.findtext("VALUE") or "").strip()
            attr_out.write(f"{s}\t{sec}\t{alias}\t{title}\t{tag}\t{val}\n")
    print(f"[ENA sample] {s} ok", flush=True)
    time.sleep(1)
attr_out.close()

# 3) GEO SOFT family files for all series + study-level ENA run list
gses = ["GSE123875", "GSE123893", "GSE123894", "GSE135230", "GSE175625", "GSE200632", "GSE293440"]
for gse in gses:
    stub = gse[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/soft/{gse}_family.soft.gz"
    data = fetch(url)
    open(f"{RAW}/{gse}_family.soft_{DATE}.gz", "wb").write(data)
    print(f"[GEO SOFT] {gse} ok ({len(data)} bytes)", flush=True)
    time.sleep(1)

studies = sorted({r[idx["secondary_study_accession"]] for r in rows})
for st in studies:
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={st}&result=read_run&fields={FIELDS}&format=tsv"
    txt = fetch(url).decode()
    open(f"{RAW}/ena_study_runs_{st}_{DATE}.tsv", "w").write(txt)
    print(f"[ENA study] {st} ok ({len(txt.splitlines())-1} runs)", flush=True)
    time.sleep(1)
print("A1 DONE", flush=True)
