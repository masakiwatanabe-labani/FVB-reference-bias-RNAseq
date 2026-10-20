#!/usr/bin/env python3
"""Phase B-4: search ENA (primary) for mouse RNA-seq studies whose SAMPLE strain attribute contains both
FVB (FVB/NJ, FVB/N or FVB) and C57BL/6J (or C57BL/6) animals. Matching is on the BioSample strain
attribute, never on study titles or background-strain strings. Output: revision_v2/audit/candidate_datasets.tsv
(+ raw query results in revision_v2/audit/raw/). GEO eutils is tried as a secondary cross-check.
"""
import urllib.request, urllib.parse, time, datetime, collections, os, re
PROJ="/usr/local/jupyter/FVB_bias"; RAW=f"{PROJ}/revision_v2/audit/raw"; AUD=f"{PROJ}/revision_v2/audit"
DATE=datetime.date.today().isoformat()
KNOWN={"SRP173496","SRP173500","SRP173501","SRP217051","SRP321579","SRP369603","SRP573960"}
FVB=["FVB/NJ","FVB/N","FVB","FVB/NCrl","FVB/NJ (JAX)","FVB/NHsd"]; B6=["C57BL/6J","C57BL/6","C57BL/6NJ","C57BL/6N","C57BL6J","C57BL6","B6"]
FIELDS="run_accession,secondary_study_accession,study_accession,sample_accession,strain,sample_title,tissue_type,cell_type,library_source,library_layout,read_count,first_public,instrument_model"
def get(url, retries=4):
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r: return r.read().decode()
        except Exception as e:
            print("  retry", i+1, e, flush=True); time.sleep(3*(i+1))
    return ""
def ena(query, fields=FIELDS, result="read_run"):
    q=urllib.parse.quote(query)
    return get(f"https://www.ebi.ac.uk/ena/portal/api/search?result={result}&query={q}&fields={fields}&limit=0&format=tsv")
def parse(t):
    L=t.splitlines()
    if not L: return []
    h=L[0].split("\t"); return [dict(zip(h,l.split("\t"))) for l in L[1:]]

# 1) all mouse RNA-seq runs with an FVB strain attribute
fvb_runs=[]
for s in FVB:
    t=ena(f'tax_eq(10090) AND library_strategy="RNA-Seq" AND strain="{s}"')
    open(f"{RAW}/ena_search_strain_{s.replace('/','_').replace(' ','_')}_{DATE}.tsv","w").write(t)
    rows=parse(t); print(f'strain="{s}": {len(rows)} runs', flush=True); fvb_runs+=rows; time.sleep(1)
studies=collections.defaultdict(list)
for r in fvb_runs: studies[r["secondary_study_accession"]].append(r)
print("FVB studies:", len(studies), flush=True)

# 2) for each study, does it also contain C57BL/6 samples (strain attribute)?
out=open(f"{AUD}/candidate_datasets.tsv","w")
out.write("study\tbioproject\tknown\tn_runs_total\tn_fvb_runs\tfvb_strains\tn_b6_runs\tb6_strains\tn_other_runs\tother_strains\tlibrary_source\tlayout\tfirst_public\ttissue_or_cell\tfvb_titles\tb6_titles\n")
n_cand=0
for st in sorted(studies):
    import glob
    cached=sorted(glob.glob(f"{RAW}/ena_study_runs_{st}_*.tsv"))
    if cached: t=open(cached[-1]).read()
    else:
        t=ena(f'secondary_study_accession="{st}"'); time.sleep(0.5)
        if t and t.startswith("run_accession"): open(f"{RAW}/ena_study_runs_{st}_{DATE}.tsv","w").write(t)
    rows=parse(t)
    if not rows or "strain" not in rows[0]: print("  skip (no strain field / empty):", st, flush=True); continue
    strains=collections.Counter(r.get("strain","") for r in rows)
    fvb=[r for r in rows if r.get("strain","") in FVB]; b6=[r for r in rows if r.get("strain","") in B6]
    other=[r for r in rows if r.get("strain","") not in FVB+B6]
    if not b6: continue
    n_cand+=1
    tis=collections.Counter((r["tissue_type"] or r["cell_type"]) for r in rows)
    out.write("\t".join([st, rows[0]["study_accession"], "Y" if st in KNOWN else "N", str(len(rows)), str(len(fvb)),
        ";".join(f"{k}:{v}" for k,v in collections.Counter(r['strain'] for r in fvb).items()), str(len(b6)),
        ";".join(f"{k}:{v}" for k,v in collections.Counter(r['strain'] for r in b6).items()), str(len(other)),
        ";".join(f"{k}:{v}" for k,v in collections.Counter(r['strain'] for r in other).items())[:200],
        ";".join(sorted(set(r["library_source"] for r in rows))), ";".join(sorted(set(r["library_layout"] for r in rows))),
        min(r["first_public"] for r in rows), ";".join(f"{k}:{v}" for k,v in tis.most_common(6))[:300],
        " | ".join(sorted(set(r["sample_title"] for r in fvb))[:6])[:400], " | ".join(sorted(set(r["sample_title"] for r in b6))[:6])[:400]])+"\n")
    out.flush(); print(f"candidate {st}: fvb={len(fvb)} b6={len(b6)} total={len(rows)}", flush=True)
out.close(); print("ENA candidates:", n_cand, flush=True)

# 3) GEO eutils cross-check (may be rate-limited / blocked)
try:
    term=urllib.parse.quote('("FVB/NJ"[All Fields] OR "FVB/N"[All Fields]) AND ("C57BL/6J"[All Fields] OR "C57BL/6"[All Fields]) AND "Mus musculus"[Organism] AND "expression profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]')
    t=get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={term}&retmax=500&retmode=json")
    open(f"{RAW}/geo_esearch_fvb_b6_{DATE}.json","w").write(t)
    import json; ids=json.loads(t)["esearchresult"]["idlist"]; print("GEO esearch GDS ids:", len(ids), flush=True)
    for i in range(0,len(ids),100):
        tt=get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={','.join(ids[i:i+100])}&retmode=json")
        open(f"{RAW}/geo_esummary_fvb_b6_{i}_{DATE}.json","w").write(tt); time.sleep(1)
except Exception as e:
    print("GEO eutils failed:", e, flush=True)
print("B4 DONE", flush=True)
