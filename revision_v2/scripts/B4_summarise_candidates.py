#!/usr/bin/env python3
"""Phase B-4 (offline step): from the cached ENA per-study run lists (revision_v2/audit/raw/ena_study_runs_*.tsv)
of every mouse RNA-seq study that has at least one run with an FVB strain attribute, find studies that also
contain C57BL/6(J) animals by the sample strain attribute (normalised, case-insensitive). Congenic / transgenic
background strings (e.g. "B6.FVB-Tg", "FVB/N-Tg(...)") are flagged and never counted as plain FVB or B6.
Writes revision_v2/audit/candidate_datasets.tsv and candidate_datasets.md."""
import glob, csv, re, collections, datetime
PROJ="/usr/local/jupyter/FVB_bias"; RAW=f"{PROJ}/revision_v2/audit/raw"; AUD=f"{PROJ}/revision_v2/audit"
KNOWN={"SRP173496":"GSE123875 (used)","SRP173500":"GSE123893 (used)","SRP173501":"GSE123894 (used)","SRP217051":"GSE135230 (old, excluded)",
       "SRP321579":"GSE175625 (old, PyVT)","SRP369603":"GSE200632 (F1)","SRP573960":"GSE293440 (F1)"}
def norm(s): return re.sub(r"\s+","",s.lower())
def cls(strain):
    n=norm(strain)
    if not n: return "none"
    if "tg(" in n or re.match(r"^(b6|c57bl/6j?)\.", n) or re.match(r"^fvb(/n|/nj)?\.", n) or "x" in n.replace("x-",""): return "background/cross"
    if n in ("fvb","fvb/n","fvb/nj","fvb/ncrl","fvb/nhsd","fvb/njrj","fvb/ntac"): return "FVB"
    if n in ("c57bl/6j","c57bl/6","c57bl6j","c57bl6","b6","c57bl/6jrj","c57bl/6jjcl","c57bl/6n","c57bl/6nj","c57bl/6ncrl","c57bl/6ntac"): return "B6"
    return "other"
fvb_studies=set()
for f in glob.glob(f"{RAW}/ena_search_strain_FVB*_*.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"): fvb_studies.add(r["secondary_study_accession"])
rows_out=[]
for st in sorted(fvb_studies):
    files=sorted(glob.glob(f"{RAW}/ena_study_runs_{st}_*.tsv"))
    if not files: rows_out.append((st,"NOT FETCHED")); continue
    rows=list(csv.DictReader(open(files[-1]), delimiter="\t"))
    if not rows or "strain" not in rows[0]: rows_out.append((st,"no strain field")); continue
    c=collections.Counter(cls(r.get("strain") or "") for r in rows)
    if c["FVB"]==0 or c["B6"]==0: continue
    fvb=[r for r in rows if cls(r.get("strain") or "")=="FVB"]; b6=[r for r in rows if cls(r.get("strain") or "")=="B6"]
    tis=collections.Counter(((r.get("tissue_type") or r.get("cell_type") or "")) for r in fvb+b6)
    rows_out.append((st, dict(bioproject=rows[0]["study_accession"], known=KNOWN.get(st,""), n_total=len(rows), n_fvb=len(fvb), n_b6=len(b6),
        fvb_strains=";".join(sorted(set((r.get("strain") or "") for r in fvb))), b6_strains=";".join(sorted(set((r.get("strain") or "") for r in b6))),
        n_background=c["background/cross"], n_other=c["other"], source=";".join(sorted(set((r.get("library_source") or "") for r in fvb+b6))),
        layout=";".join(sorted(set((r.get("library_layout") or "") for r in fvb+b6))), first_public=min((r.get("first_public") or "") for r in rows),
        tissues=";".join(f"{k}:{v}" for k,v in tis.most_common(5)), fvb_titles=" | ".join(sorted(set((r.get("sample_title") or "") for r in fvb))[:5])[:300],
        b6_titles=" | ".join(sorted(set((r.get("sample_title") or "") for r in b6))[:5])[:300])))
cands=[(s,d) for s,d in rows_out if isinstance(d,dict)]
with open(f"{AUD}/candidate_datasets.tsv","w") as o:
    keys=["bioproject","known","n_total","n_fvb","n_b6","fvb_strains","b6_strains","n_background","n_other","source","layout","first_public","tissues","fvb_titles","b6_titles"]
    o.write("study\t"+"\t".join(keys)+"\n")
    for s,d in sorted(cands, key=lambda x:(-min(x[1]["n_fvb"],x[1]["n_b6"]))): o.write(s+"\t"+"\t".join(str(d[k]) for k in keys)+"\n")
with open(f"{AUD}/candidate_datasets.md","w") as o:
    o.write(f"# Phase B-4 — search for additional FVB vs C57BL/6J RNA-seq datasets ({datetime.date.today()})\n\n")
    o.write("Search: ENA portal API, `result=read_run`, `tax_eq(10090) AND library_strategy=\"RNA-Seq\" AND strain=\"<s>\"` for s in FVB/NJ, FVB/N, FVB, FVB/NCrl, FVB/NHsd (raw results `revision_v2/audit/raw/ena_search_strain_*`). "
            f"{len(fvb_studies)} studies had ≥1 FVB-attributed run; each study's full run list was fetched (`ena_study_runs_<study>_*.tsv`) and the sample **strain attribute** of every run was classified (plain FVB / plain B6 / congenic-transgenic background / other). "
            "A study is a candidate only if it contains both plain-FVB and plain-B6 animals; strings such as `B6.FVB-Tg(...)` or `FVB/N-Tg(Thy1-cre)` are classed as background and never count (rule B-4). GEO eutils was queried as a cross-check (see log).\n\n")
    o.write(f"Not fetched / no strain field: {sum(1 for s,d in rows_out if not isinstance(d,dict))} studies (listed at the end).\n\n")
    o.write(f"## Candidates with both strains as sample attributes: {len(cands)}\n\n| study | BioProject | known | runs FVB / B6 / total | FVB strain strings | B6 strain strings | source | layout | first public | tissue/cell (top) |\n|---|---|---|---|---|---|---|---|---|---|\n")
    for s,d in sorted(cands, key=lambda x:(-min(x[1]["n_fvb"],x[1]["n_b6"]))):
        o.write(f"| {s} | {d['bioproject']} | {d['known']} | {d['n_fvb']} / {d['n_b6']} / {d['n_total']} | {d['fvb_strains']} | {d['b6_strains']} | {d['source']} | {d['layout']} | {d['first_public']} | {d['tissues'][:120]} |\n")
    o.write("\n## Sample titles (for manual screening of tissue / genotype / treatment)\n\n")
    for s,d in sorted(cands, key=lambda x:(-min(x[1]["n_fvb"],x[1]["n_b6"]))):
        o.write(f"- **{s}**: FVB: {d['fvb_titles']} — B6: {d['b6_titles']}\n")
    o.write("\n## Not evaluated\n\n" + "\n".join(f"- {s}: {d}" for s,d in rows_out if not isinstance(d,dict)) + "\n")
print("FVB studies:", len(fvb_studies), "candidates:", len(cands))
