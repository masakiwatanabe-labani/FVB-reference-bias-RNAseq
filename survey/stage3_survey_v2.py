#!/usr/bin/env python3
"""
Stage 3 Survey v2 — honest classification with evidence
========================================================
Categories (mutually exclusive, priority order):
  strain_aware : mentions WASP / g2gtools / FVB genome / personalized ref / N-masking
  stated_B6    : explicitly mentions GRCm38/39 / mm10/mm39 / C57BL/6 / Bowtie+mm / STAR+GRCm
  not_reported : data_processing fetched but no reference genome keyword found
  fetch_failed : GEO returned no data_processing text (network error or field absent)

Each row stores:
  evidence_text   : the actual data_processing sentences that triggered classification
  evidence_source : "GEO_sample_dp" | "GEO_series_dp" | "none"

Fig 1 should show stated_B6 / not_reported / fetch_failed / strain_aware separately.
NOT silently merging not_reported + fetch_failed into B6.
"""

import sys, os, time, re, csv
from urllib.request import urlopen, Request
from urllib.error import URLError

# Usage: stage3_survey_v2.py [INPUT_TSV] [OUTPUT_DIR]
#   INPUT_TSV : list of GEO series with a column gse_acc (and optionally title); default survey_classification.tsv
#   OUTPUT_DIR: default ./survey_output (survey_of_practice_v2.tsv, summary_stats_v2.txt)
INPUT_TSV  = sys.argv[1] if len(sys.argv) > 1 else "survey_classification.tsv"
OUT_DIR    = sys.argv[2] if len(sys.argv) > 2 else "survey_output"
OUTPUT_TSV = os.path.join(OUT_DIR, "survey_of_practice_v2.tsv")
DELAY      = 1.2   # NCBI polite rate: <3 req/sec without API key

if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

# ── keyword lists ────────────────────────────────────────────────────────────
STRAIN_AWARE_KW = [
    "wasp", "g2gtools", "personali", "strain-aware", "strain aware",
    "fvb genome", "fvb/n genome", "fvb/nj genome",
    "allele-specific alignment", "snp-aware", "n-masking", "n masking",
    "diploid genome", "individualized genome",
]
B6_KW = [
    "grcm39", "grcm38", "mm10", "mm39", "mm9",
    "c57bl/6", "c57bl6", "b6 reference", "b6 genome",
    "mus_musculus.grcm", "ucsc mm", "ensembl grch",  # common alignment pipeline strings
    "genome_mouse", "ncbi37", "ncbi38",
    # aligner + implicit B6 (common in STAR/HISAT2 commands with GRCm paths)
    "/grcm3", "/mm10/", "/mm39/",
]

def fetch_soft(gse: str, targ: str = "gsm", view: str = "quick") -> str:
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
           f"?acc={gse}&targ={targ}&form=text&view={view}")
    try:
        req = Request(url, headers={"User-Agent": "FVB_bias_survey/2.0 (research)"})
        with urlopen(req, timeout=30) as r:
            return r.read().decode("utf-8", errors="replace")
    except URLError as e:
        print(f"    WARN fetch {gse} targ={targ}: {e}", file=sys.stderr)
        return ""

def extract_dp_lines(text: str, prefix: str) -> list[str]:
    """Extract !Sample_data_processing or !Series_data_processing lines."""
    lines = []
    for line in text.splitlines():
        if line.startswith(prefix):
            val = line.split("=", 1)[-1].strip()
            if val:
                lines.append(val)
    return lines

def classify(dp_text: str) -> tuple[str, str]:
    """
    Returns (category, matched_keyword).
    Priority: strain_aware > stated_B6 > not_reported
    """
    t = dp_text.lower()
    for kw in STRAIN_AWARE_KW:
        if kw in t:
            return "strain_aware", kw
    for kw in B6_KW:
        if kw in t:
            return "stated_B6", kw
    return "not_reported", ""

def survey_one(gse: str) -> dict:
    result = {
        "gse_acc": gse,
        "category": "fetch_failed",
        "evidence_text": "",
        "evidence_source": "none",
        "matched_keyword": "",
    }

    # 1) Try targ=gsm (sample-level data_processing — most reliable)
    gsm_text = fetch_soft(gse, targ="gsm", view="quick")
    time.sleep(DELAY)
    dp_lines = extract_dp_lines(gsm_text, "!Sample_data_processing")

    if dp_lines:
        dp_combined = " | ".join(dict.fromkeys(dp_lines))  # deduplicate
        cat, kw = classify(dp_combined)
        result.update({
            "category": cat,
            "evidence_text": dp_combined[:600],
            "evidence_source": "GEO_sample_dp",
            "matched_keyword": kw,
        })
        return result

    # 2) Fallback: targ=self (series-level)
    series_text = fetch_soft(gse, targ="self", view="quick")
    time.sleep(DELAY)
    dp_lines2 = extract_dp_lines(series_text, "!Series_data_processing")

    if dp_lines2:
        dp_combined2 = " | ".join(dict.fromkeys(dp_lines2))
        cat, kw = classify(dp_combined2)
        result.update({
            "category": cat,
            "evidence_text": dp_combined2[:600],
            "evidence_source": "GEO_series_dp",
            "matched_keyword": kw,
        })
        return result

    # 3) No data_processing found in GEO at all
    result["category"] = "not_reported"
    result["evidence_source"] = "GEO_fetched_no_dp"
    return result

def main():
    studies = []
    with open(INPUT_TSV, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            studies.append(row)
    print(f"Loaded {len(studies)} studies")

    fieldnames = [
        "gse_acc", "title", "category",
        "evidence_source", "matched_keyword", "evidence_text",
        "validation_url",
    ]

    # Resume support: skip already-done
    done = set()
    rows_out = []
    if os.path.exists(OUTPUT_TSV):
        with open(OUTPUT_TSV, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                done.add(row["gse_acc"])
                rows_out.append(row)
        print(f"Resuming: {len(done)} already done")

    with open(OUTPUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows_out)

        for i, study in enumerate(studies):
            gse   = study["gse_acc"].strip()
            title = study.get("title", "").strip()
            if gse in done:
                continue

            print(f"[{i+1}/{len(studies)}] {gse} ...", flush=True)
            res = survey_one(gse)
            res["title"] = title
            res["validation_url"] = (
                f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}")

            print(f"  → {res['category']} | src={res['evidence_source']} "
                  f"| kw={res['matched_keyword']!r} "
                  f"| dp={res['evidence_text'][:60]!r}")

            writer.writerow(res)
            f.flush()

    # ── Summary ────────────────────────────────────────────────────────────────
    all_rows = []
    with open(OUTPUT_TSV, newline="") as f:
        all_rows = list(csv.DictReader(f, delimiter="\t"))

    total = len(all_rows)
    from collections import Counter
    cats = Counter(r["category"] for r in all_rows)

    print("\n=== SURVEY SUMMARY v2 ===")
    print(f"Total: {total}")
    for cat in ["stated_B6","not_reported","strain_aware","fetch_failed"]:
        n = cats.get(cat, 0)
        print(f"  {cat:20s}: {n:3d} ({100*n/total:.1f}%)")
    print(f"\n  Evidence coverage: "
          f"{sum(1 for r in all_rows if r['evidence_text'])}/{total} have dp text")
    print(f"  fetch_failed (no dp found at all): {cats.get('fetch_failed',0)}")

    summary = os.path.join(OUT_DIR, "summary_stats_v2.txt")
    with open(summary, "w") as f:
        f.write(f"Total studies: {total}\n")
        for cat in ["stated_B6","not_reported","strain_aware","fetch_failed"]:
            n = cats.get(cat, 0)
            f.write(f"{cat}: {n} ({100*n/total:.1f}%)\n")
        f.write(f"evidence_coverage: "
                f"{sum(1 for r in all_rows if r['evidence_text'])}/{total}\n")
    print(f"Summary → {summary}")
    print(f"TSV     → {OUTPUT_TSV}")

if __name__ == "__main__":
    main()
