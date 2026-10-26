#!/usr/bin/env python3
"""Phase W-1 re-screen, step 1: fetch and cache the raw GEO/SRA records of every candidate series.
For each GSE in the input list (column gse_acc):
  raw/<GSE>.series.soft.gz   GEO SOFT, series record (targ=self, view=brief): !Series_type, summary, overall design, relations
  raw/<GSE>.samples.soft.gz  GEO SOFT, all sample records (targ=gsm, view=brief): characteristics, library_strategy/source,
                             data_processing, SRA/BioSample relations
  raw/<GSE>.sra_runinfo.csv.gz  SRA run table (E-utilities efetch rettype=runinfo) of the SRX experiments linked from the samples:
                             LibraryStrategy, LibrarySource, SampleName (GSM), BioSample
The submitters' contact lines (!Series_contact_*, !Sample_contact_*: name, e-mail, address, phone) are not needed for the re-screen
and are removed before caching (also from cache files written before this filter existed). Existing cache files are not re-fetched.
Usage: rescreen_fetch.py INPUT_TSV RAW_DIR"""
import sys, os, time, gzip, re, csv, threading
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from urllib.error import URLError, HTTPError

INPUT_TSV, RAW = sys.argv[1], sys.argv[2]
DELAY = 0.4
UA = {"User-Agent": "FVB_bias_survey/2.0 (research)"}
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

_lock = threading.Lock(); _last = [0.0]
def get(url, data=None, tries=4):
    for k in range(tries):
        if url.startswith(EU):   # E-utilities: at most ~2.5 requests/s over all threads (no API key)
            with _lock:
                time.sleep(max(0, _last[0] + DELAY - time.time())); _last[0] = time.time()
        try:
            req = Request(url, data=urlencode(data).encode() if data else None, headers=UA)
            with urlopen(req, timeout=120) as r:
                out = r.read().decode("utf-8", errors="replace")
            return out
        except (URLError, HTTPError, TimeoutError) as e:
            print(f"  retry {k+1}: {e} [{url[:90]}]", file=sys.stderr); time.sleep(3 * (k + 1))
    raise RuntimeError(f"failed: {url}")

def cache(path, fn):
    if os.path.exists(path): return gzip.open(path, "rt", encoding="utf-8").read()
    txt = fn()
    with gzip.open(path + ".part", "wt", encoding="utf-8") as o: o.write(txt)
    os.replace(path + ".part", path)   # atomic: an interrupted fetch leaves no truncated cache file
    return txt

CONTACT = re.compile(r"^!(Series|Sample|Platform)_contact_[^\n]*\n?", re.M)
def soft(gse, targ):
    return CONTACT.sub("", get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ={targ}&form=text&view=brief"))

def runinfo(srx):
    if not srx: return ""
    rows, header = [], None
    for i in range(0, len(srx), 150):
        chunk = srx[i:i + 150]
        es = get(EU + "esearch.fcgi", {"db": "sra", "term": " OR ".join(f"{x}[Accession]" for x in chunk),
                                       "retmax": 10000, "usehistory": "y"})
        we = re.search(r"<WebEnv>(\S+?)</WebEnv>", es); qk = re.search(r"<QueryKey>(\d+)</QueryKey>", es)
        if not we: continue
        txt = get(EU + "efetch.fcgi", {"db": "sra", "query_key": qk.group(1), "WebEnv": we.group(1), "rettype": "runinfo", "retmode": "text"})
        lines = [l for l in txt.splitlines() if l.strip()]
        if not lines: continue
        if header is None: header = lines[0]
        rows += [l for l in lines[1:] if l != lines[0]]
    return "\n".join([header] + rows) + "\n" if header else ""

if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    gses = [r["gse_acc"] for r in csv.DictReader(open(INPUT_TSV), delimiter="\t")]
    LOG = open(f"{RAW}/fetch.log", "a"); T0 = time.time(); done = [0]; lk = threading.Lock()
    LOG.write(f"# start {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: {len(gses)} series, 4 workers\n"); LOG.flush()
    def one(ig):
        i, g = ig
        s = cache(f"{RAW}/{g}.series.soft.gz", lambda: soft(g, "self"))
        m = cache(f"{RAW}/{g}.samples.soft.gz", lambda: soft(g, "gsm"))
        srx = sorted(set(re.findall(r"^!Sample_relation = SRA: \S*?term=(SRX\d+)", m, flags=re.M)))
        ri = cache(f"{RAW}/{g}.sra_runinfo.csv.gz", lambda: runinfo(srx))
        msg = f"[{i}/{len(gses)}] {g}: samples={m.count('^SAMPLE')} srx={len(srx)} runs={max(0, ri.count(chr(10)) - 1)}"
        with lk:
            done[0] += 1; el = time.time() - T0
            LOG.write(f"{time.strftime('%H:%M:%S')} {msg} elapsed={el:.0f}s\n")
            if done[0] % 100 == 0: LOG.write(f"## {done[0]} series done, {el/60:.1f} min elapsed ({el/done[0]*100/60:.1f} min per 100)\n")
            LOG.flush()
        print(msg, flush=True)
    with ThreadPoolExecutor(4) as ex: list(ex.map(one, enumerate(gses, 1)))   # 4 series at a time (GEO acc.cgi is slow)
    for g in gses:   # contact lines out of every cached SOFT record (idempotent)
        for kind in ("series", "samples"):
            f = f"{RAW}/{g}.{kind}.soft.gz"; t = gzip.open(f, "rt", encoding="utf-8").read(); u = CONTACT.sub("", t)
            if u != t:
                with gzip.open(f + ".part", "wt", encoding="utf-8") as o: o.write(u)
                os.replace(f + ".part", f)
    if not os.path.exists(f"{RAW}/FETCHED.txt"):
      with open(f"{RAW}/FETCHED.txt", "w") as o:
        o.write(f"fetched {time.strftime('%Y-%m-%d %H:%M %Z')} from GEO (acc.cgi, form=text, view=brief) and SRA (E-utilities runinfo); {len(gses)} series\n")
