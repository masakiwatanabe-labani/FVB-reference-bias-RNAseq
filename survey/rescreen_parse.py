"""Phase W-1 re-screen: parsers for the cached GEO SOFT and SRA run-table records (raw/, written by rescreen_fetch.py)."""
import gzip, re, csv, io, collections

def _read(path):
    try: return gzip.open(path, "rt", encoding="utf-8").read()
    except FileNotFoundError: return ""

def series(raw, gse):
    """series record: title, types, summary, overall_design, superseries_of, subseries_of"""
    t = _read(f"{raw}/{gse}.series.soft.gz"); d = collections.defaultdict(list)
    for l in t.splitlines():
        if l.startswith("!Series_"):
            k, _, v = l.partition(" = "); d[k[8:]].append(v.strip())
    rel = d.get("relation", [])
    return dict(title=(d.get("title") or [""])[0], types=d.get("type", []), summary=" ".join(d.get("summary", [])),
                overall_design=" ".join(d.get("overall_design", [])),
                superseries_of=[x.split(":", 1)[1].strip() for x in rel if x.startswith("SuperSeries of:")],
                subseries_of=[x.split(":", 1)[1].strip() for x in rel if x.startswith("SubSeries of:")],
                sample_ids=d.get("sample_id", []))

def samples(raw, gse):
    """list of sample dicts: gsm, title, source, organism, characteristics (list of 'key: value'), strategy, source_lib, molecule,
    data_processing (list), srx, type; library strategy/source from the SRA run table where the SRX is in it, else the GEO field"""
    t = _read(f"{raw}/{gse}.samples.soft.gz"); out = []; cur = None
    for l in t.splitlines():
        if l.startswith("^SAMPLE"):
            cur = dict(gsm=l.split("=", 1)[1].strip(), title="", source="", organism=[], characteristics=[], strategy="", source_lib="",
                       molecule="", data_processing=[], srx="", type="", description=[], extract=[], protocol=[]); out.append(cur); continue
        if cur is None or not l.startswith("!Sample_"): continue
        k, _, v = l.partition(" = "); k = k[8:]; v = v.strip()
        if k == "title": cur["title"] = v
        elif k.startswith("source_name"): cur["source"] = (cur["source"] + " | " + v).strip(" |")
        elif k.startswith("organism"): cur["organism"].append(v)
        elif k.startswith("characteristics"): cur["characteristics"].append(v)
        elif k == "library_strategy": cur["strategy"] = v
        elif k == "library_source": cur["source_lib"] = v
        elif k.startswith("molecule"): cur["molecule"] = v
        elif k == "data_processing": cur["data_processing"].append(v)
        elif k == "type": cur["type"] = v
        elif k == "description": cur["description"].append(v)
        elif k.startswith("extract_protocol"): cur["extract"].append(v)
        elif k.startswith(("growth_protocol", "treatment_protocol")): cur["protocol"].append(v)
        elif k == "relation" and v.startswith("SRA:"):
            m = re.search(r"(SRX\d+)", v); cur["srx"] = m.group(1) if m else ""
    ri = _read(f"{raw}/{gse}.sra_runinfo.csv.gz"); sra = {}
    if ri.strip():
        for r in csv.DictReader(io.StringIO(ri)): sra[r["Experiment"]] = r
    for s in out:
        r = sra.get(s["srx"])
        s["strategy_geo"], s["source_lib_geo"] = s["strategy"], s["source_lib"]
        if r:
            s["strategy"], s["source_lib"], s["strategy_from"] = r["LibraryStrategy"], r["LibrarySource"], "SRA"
        else:
            s["strategy_from"] = "GEO" if s["strategy"] else "none"
    return out
