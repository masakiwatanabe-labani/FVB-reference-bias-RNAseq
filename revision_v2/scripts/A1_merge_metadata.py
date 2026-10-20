#!/usr/bin/env python3
"""Phase A-1: merge ENA run report + ENA/BioSample attributes + GEO SOFT into one per-SRR table
(revision_v2/audit/old_runs_metadata.tsv) and a per-series summary (series_structure.tsv).
"""
import gzip, glob, os, re, collections
PROJ = "/usr/local/jupyter/FVB_bias"; AUD = f"{PROJ}/revision_v2/audit"; RAW = f"{AUD}/raw"

# --- old group assignment
old = {}
for l in open(f"{AUD}/old_group_assignment.tsv").read().splitlines()[1:]:
    srr, gse, grp, src = l.split("\t"); old[srr] = (gse, grp)

# --- ENA runs
rows = [l.split("\t") for l in open(f"{AUD}/old_runs_ena.tsv").read().splitlines()]
h = {k: i for i, k in enumerate(rows[0])}
ena = {r[h["run_accession"]]: r for r in rows[1:]}

# --- BioSample attributes (ENA XML)
attrs = collections.defaultdict(dict); samp_title = {}; samp_alias = {}; samp_sec = {}
for l in open(f"{AUD}/old_samples_attributes.tsv").read().splitlines()[1:]:
    s, sec, alias, title, tag, val = (l.split("\t") + [""] * 6)[:6]
    attrs[s][tag.lower()] = val; samp_title[s] = title; samp_alias[s] = alias; samp_sec[s] = sec

# --- GEO SOFT
geo_samples = {}; series_info = {}
for p in sorted(glob.glob(f"{RAW}/GSE*_family.soft_*.gz")):
    cur = None; ser = None
    with gzip.open(p, "rt", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("^SERIES"):
                ser = line.split("=")[1].strip(); series_info[ser] = collections.defaultdict(list); cur = None
            elif line.startswith("^SAMPLE"):
                cur = line.split("=")[1].strip(); geo_samples[cur] = collections.defaultdict(list)
            elif line.startswith("!Series_") and cur is None and ser:
                k, _, v = line[1:].partition(" = "); series_info[ser][k].append(v)
            elif line.startswith("!Sample_") and cur:
                k, _, v = line[1:].partition(" = "); geo_samples[cur][k].append(v)

# index GEO samples by BioSample and SRX
by_biosample = {}; by_srx = {}
for gsm, d in geo_samples.items():
    for rel in d.get("Sample_relation", []):
        m = re.search(r"(SAMN\d+)", rel)
        if m: by_biosample[m.group(1)] = gsm
        m = re.search(r"(SRX\d+)", rel)
        if m: by_srx[m.group(1)] = gsm

def pick(d, *keys):
    for k in keys:
        if k in d and d[k]: return d[k]
    return ""

def geo_char(d, key):
    for c in d.get("Sample_characteristics_ch1", []):
        k, _, v = c.partition(": ")
        if k.strip().lower() == key: return v.strip()
    return ""

# technical runs per BioSample
runs_per_sample = collections.Counter(r[h["sample_accession"]] for r in rows[1:])

cols = ["SRR","GSE","old_group","GSM","BioSample","SRS","SRX","SRP","sample_title_ena","geo_title","geo_source_name",
        "geo_strain","geo_genotype","geo_tissue","geo_treatment","geo_sex","geo_age","geo_characteristics_all",
        "biosample_strain","biosample_attrs_all","library_layout","instrument","read_count","n_technical_runs","geo_series_ids"]
out = open(f"{AUD}/old_runs_metadata.tsv", "w"); out.write("\t".join(cols) + "\n")
for srr in open(f"{AUD}/old_srr_list.txt").read().split():
    r = ena[srr]; s = r[h["sample_accession"]]; srx = r[h["experiment_accession"]]
    gsm = by_biosample.get(s) or by_srx.get(srx) or ""
    g = geo_samples.get(gsm, {})
    a = attrs.get(s, {})
    row = [srr, old[srr][0], old[srr][1], gsm, s, r[h["secondary_sample_accession"]], srx, r[h["secondary_study_accession"]],
           r[h["sample_title"]], " ".join(g.get("Sample_title", [])), " ".join(g.get("Sample_source_name_ch1", [])),
           geo_char(g, "strain") or geo_char(g, "background strain") or geo_char(g, "strain background"), geo_char(g, "genotype") or geo_char(g, "genotype/variation"),
           geo_char(g, "tissue") or geo_char(g, "cell type"), geo_char(g, "treatment") or geo_char(g, "diet"),
           geo_char(g, "sex") or geo_char(g, "gender"), geo_char(g, "age"),
           " | ".join(g.get("Sample_characteristics_ch1", [])),
           pick(a, "strain", "strain background"), " | ".join(f"{k}={v}" for k, v in sorted(a.items()) if k not in ("ena-first-public","ena-last-update")),
           r[h["library_layout"]], r[h["instrument_model"]], r[h["read_count"]], str(runs_per_sample[s]),
           ",".join(g.get("Sample_series_id", []))]
    out.write("\t".join(row) + "\n")
out.close()

with open(f"{AUD}/series_structure.tsv", "w") as o:
    o.write("GSE\ttitle\trelations\tn_samples_in_soft\tsample_ids\n")
    for ser, d in series_info.items():
        o.write(f"{ser}\t{' '.join(d.get('Series_title', []))}\t{' ; '.join(d.get('Series_relation', []))}\t{len(d.get('Series_sample_id', []))}\t{','.join(d.get('Series_sample_id', []))}\n")
print("merged", sum(1 for _ in open(f"{AUD}/old_runs_metadata.tsv")) - 1, "runs")
