#!/usr/bin/env python3
"""Phase X-1: candidate GEO series for the survey from a broad E-utilities search of db=gds.
Queries (recorded with date and hit count in OUT_DIR/search_log.tsv):
  Q1  (FVB[All Fields]) AND "Mus musculus"[Organism] AND gse[Entry Type] AND "Expression profiling by high throughput sequencing"[DataSet Type]
  Q1o same with "Other"[DataSet Type];  Q1n same with "Non-coding RNA profiling by high throughput sequencing"[DataSet Type]
  Q2  (FVB[All Fields]) AND "Mus musculus"[Organism] AND gsm[Entry Type] -> GSE(s) of each GSM (esummary field gse); a GSE is kept
      when its series type includes "high throughput sequencing"
Union of the GSE accessions -> OUT_DIR/candidates_broad.tsv (gse_acc, source_query, series_type, title, release_date, n_samples).
Optional: --add-variant TERM (Phase X, authors' decision: FVBN) runs the same queries with TERM in place of FVB and appends the
accessions not yet in OUT_DIR/candidates_broad.tsv with source_query = TERM (the log gets the counts).
Optional: --variants tests the strain-name variants (FVB/N, FVB/NJ, ...) in place of FVB in Q1 and Q2 and reports added accessions.
Usage: broad_search.py OUT_DIR [--variants]"""
import sys, os, json, time, csv, urllib.request, urllib.parse
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"; DELAY = 0.4; UA = {"User-Agent": "FVB_bias_survey/2.0 (research)"}
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
def call(fn, params, tries=6):
    for k in range(tries):
        try:
            time.sleep(DELAY)
            req = urllib.request.Request(EU + fn, data=urllib.parse.urlencode(params).encode(), headers=UA)
            return json.loads(urllib.request.urlopen(req, timeout=120).read())
        except Exception as e:
            print(f"  retry {fn} {k+1}: {e}", file=sys.stderr); time.sleep(3 * (k + 1))
    raise RuntimeError(fn)
def search(term):
    r = call("esearch.fcgi", {"db": "gds", "term": term, "retmax": 0, "usehistory": "y", "retmode": "json"})["esearchresult"]
    return int(r["count"]), r["webenv"], r["querykey"]
def summaries(count, we, qk, step=500):
    out = []
    for s in range(0, count, step):
        r = call("esummary.fcgi", {"db": "gds", "query_key": qk, "WebEnv": we, "retstart": s, "retmax": step, "retmode": "json"})["result"]
        out += [r[u] for u in r.get("uids", [])]
    return out
def gse_summaries(nums):
    out = {}
    nums = sorted(set(nums))
    for i in range(0, len(nums), 400):
        r = call("esummary.fcgi", {"db": "gds", "id": ",".join(str(200000000 + int(n)) for n in nums[i:i + 400]), "retmode": "json"})["result"]
        for u in r.get("uids", []):
            x = r[u]
            if x.get("entrytype") == "GSE": out[x["accession"]] = x
    return out
BASE = '({s}[All Fields]) AND "Mus musculus"[Organism] AND {e}[Entry Type]'
TYPES = {"Q1": "Expression profiling by high throughput sequencing", "Q1o": "Other", "Q1n": "Non-coding RNA profiling by high throughput sequencing"}
def run(strain, log):
    found = {}   # gse -> set(queries)
    info = {}
    for q, t in TYPES.items():
        term = BASE.format(s=strain, e="gse") + f' AND "{t}"[DataSet Type]'; n, we, qk = search(term)
        sm = summaries(n, we, qk); log.append((q, strain, term, n, len(sm)))
        for x in sm: found.setdefault(x["accession"], set()).add(q); info[x["accession"]] = x
    term = BASE.format(s=strain, e="gsm"); n, we, qk = search(term)
    sm = summaries(n, we, qk); gnums = set()
    for x in sm:
        for g in str(x.get("gse", "")).split(";"):
            if g.strip().isdigit(): gnums.add(int(g))
    gi = gse_summaries([g for g in gnums if f"GSE{g}" not in info])
    info.update(gi); kept = 0
    for g in gnums:
        acc = f"GSE{g}"; x = info.get(acc)
        if x and "high throughput sequencing" in x.get("gdstype", ""): found.setdefault(acc, set()).add("Q2"); kept += 1
    log.append(("Q2", strain, term, n, f"{len(sm)} GSM -> {len(gnums)} GSE -> {kept} with a high-throughput-sequencing type"))
    return found, info
if __name__ == "__main__" and "--add-variant" in sys.argv:
    term = sys.argv[sys.argv.index("--add-variant") + 1]; date = time.strftime("%Y-%m-%d %H:%M %Z"); log = []
    have = {r["gse_acc"] for r in csv.DictReader(open(f"{OUT}/candidates_broad.tsv"), delimiter="\t")}
    f2, info = run(term, log); add = sorted(set(f2) - have, key=lambda a: int(a[3:]))
    with open(f"{OUT}/candidates_broad.tsv", "a", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        for g in add: x = info[g]; w.writerow([g, term, x.get("gdstype", ""), x.get("title", ""), x.get("pdat", ""), x.get("n_samples", "")])
    with open(f"{OUT}/search_log.tsv", "a") as o:
        o.write(f"# added {date}: variant {term} (same queries, {term} in place of FVB)\n")
        for r in log: o.write("\t".join(map(str, r)) + "\n")
        o.write(f"variant_added\t{term}\t\t{len(f2)}\t{len(add)} accessions not in the FVB union appended (source_query = {term})\n")
    print("added", len(add)); sys.exit(0)
if __name__ == "__main__":
    date = time.strftime("%Y-%m-%d %H:%M %Z"); log = []
    found, info = run("FVB", log)
    with open(f"{OUT}/candidates_broad.tsv", "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["gse_acc", "source_query", "series_type", "title", "release_date", "n_samples"])
        for g in sorted(found, key=lambda a: int(a[3:])):
            x = info[g]; w.writerow([g, ";".join(sorted(found[g])), x.get("gdstype", ""), x.get("title", ""), x.get("pdat", ""), x.get("n_samples", "")])
    if "--variants" in sys.argv:
        for v in ('"FVB/N"', '"FVB/NJ"', '"FVB/NCrl"', '"FVB/NTac"', "FVBN", '"FVB-N"'):
            f2, _ = run(v, log); add = sorted(set(f2) - set(found), key=lambda a: int(a[3:]))
            log.append(("variant", v, "", len(f2), "adds " + (",".join(add) if add else "nothing")))
    with open(f"{OUT}/search_log.tsv", "w") as o:
        o.write(f"# run {date}\nquery\tstrain_term\tterm\tcount\tretrieved\n")
        for r in log: o.write("\t".join(map(str, r)) + "\n")
    print(open(f"{OUT}/search_log.tsv").read()); print("union:", len(found))
