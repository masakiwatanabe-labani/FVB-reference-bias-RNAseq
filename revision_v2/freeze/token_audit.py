#!/usr/bin/env python3
"""Phase M-1: numeric-token audit of every deliverable (standalone; reads the deliverables, numbers.tsv and the xlsx only).
Scope: revision_v2/submission/{LAR_main_text_v2,LAR_additional_file_1,response_to_reviewers_v2,cover_letter_v2}.plain.txt and
`pdftotext -layout` of revision_v2/submission/figures/*.pdf. Every numeric token is classed, in this precedence:
  STRUCT   citation/figure/table/panel labels, accessions, years, versions, parameters, reference list, title page, list ordinals,
           figure axis/colour-bar ticks, thresholds and design constants fixed in Methods (reason in `ref`);
  KEY      equals a numbers.tsv value at the printed precision (value, |value|, x100 for %, /1e6 for "million"; keys in `ref`);
  XLSX     equals, at printed precision, a numeric cell, a column min/max/mean/median, a row count or a categorical value count
           of a sheet of supplementary_tables_v2.xlsx (sheet:column in `ref`);
  UNMAPPED anything else.
Figure tokens that are not STRUCT get a `note`: the text file(s) in which the same number is printed, or "no text counterpart".
Output: revision_v2/freeze/token_audit.tsv (file, section, sentence, token, class, ref, note). Usage: token_audit.py"""
import re, csv, subprocess, collections, math, glob, os
import pandas as pd
P = "/usr/local/jupyter/FVB_bias"; SUB = f"{P}/revision_v2/submission"; OUT = f"{P}/revision_v2/freeze/token_audit.tsv"
SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")

# ---------------------------------------------------------------- reference values ------------------------------------------------
KEYS = []
for r in csv.DictReader(open(f"{P}/manuscript/revised_v2/generated/numbers.tsv"), delimiter="\t"):
    try: x = float(r["value"])
    except ValueError: continue
    if math.isnan(x) or math.isinf(x): continue
    KEYS.append((r["key"], x))
def cands(x):
    c = {abs(x), abs(x) * 100, abs(x) / 1e6}
    return c
KEY_IDX = collections.defaultdict(set)          # (decimals, string) -> keys
for k, x in KEYS:
    for c in cands(x):
        for d in range(0, 5): KEY_IDX[(d, f"{c:.{d}f}")].add(k)
KEY_SCI = collections.defaultdict(set)          # (mantissa decimals, "m e") -> keys
for k, x in KEYS:
    if x == 0: continue
    for d in range(0, 3):
        m, e = f"{abs(x):.{d}e}".split("e"); KEY_SCI[(d, m, int(e))].add(k)
KEY_ZERO_P = [k for k, x in KEYS if (x == 0 or abs(x) < 1e-15) and ("_p" in k or k.endswith("p") or "pval" in k or "p_" in k)]

XL_IDX = collections.defaultdict(set)
xl = pd.read_excel(f"{SUB}/supplementary_tables_v2.xlsx", sheet_name=None, dtype=str)
def xl_add(v, ref):
    for d in range(0, 5): XL_IDX[(d, f"{abs(v):.{d}f}")].add(ref); XL_IDX[(d, f"{abs(v) * 100:.{d}f}")].add(ref + "(%)")
for sh, df in xl.items():
    if sh == "Contents": continue
    xl_add(len(df), f"{sh}:nrows")
    for col in df.columns:
        s = df[col].replace("", None)
        num = pd.to_numeric(s, errors="coerce")
        if num.notna().sum() > 0:
            vals = num.dropna()
            for v in vals.unique(): xl_add(float(v), f"{sh}:{col}")
            for nm, v in (("min", vals.min()), ("max", vals.max()), ("mean", vals.mean()), ("median", vals.median())): xl_add(float(v), f"{sh}:{col}:{nm}")
        if num.notna().sum() < len(s.dropna()):
            for val, n in s.dropna().value_counts().items(): xl_add(float(n), f"{sh}:{col}=count({str(val)[:30]})")

# ---------------------------------------------------------------- tokenisation -------------------------------------------------------
NUM = r"(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?"
TOK = re.compile(
    r"(?P<sci>(?P<sm>[−\-+]?\d+(?:\.\d+)?)\s*×\s*10(?P<se>[⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+))"
    r"|(?P<enot>(?<![\w.])\d+(?:\.\d+)?e[-+]\d+)"
    r"|(?P<pow>(?<![\w.])10(?P<pe>[⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+))"
    r"|(?P<plain>(?<![\w.\-−+])(?P<sign>[−\-+±]?)(?<![\w.])(?P<num>" + NUM + r")(?P<pct>%?)(?![\w.]*[A-Za-z_]|\.\d|[⁰¹²³⁴⁵⁶⁷⁸⁹₀₁₂₃₄₅₆₇₈₉]))")
STRIP = re.compile(r"https?://\S+|doi\.org/\S+|\b10\.\d{4,}/\S+|\S+@\S+|`?[\w./-]+\.(?:tsv|md|py|sh|bed|txt|R|xlsx|vcf|gtf|gz|bam|docx|pdf|tiff|md5|yaml|json)\b`?")
def blank(m): return " " * len(m.group(0))

SOFTWARE = r"(?:STAR|samtools|featureCounts|Subread|edgeR|limma|g2gtools|WASP|Python|R|bcftools|bedtools|pandas|scipy|SciPy|numpy|NumPy|matplotlib|pysam|GATK|Bioconductor|fastp|Trimmomatic|sra-tools|SRA Toolkit|fasterq-dump|pandoc|version|v\.?)"
LABEL = r"(?:Lines?|Figs?\.|Figures?|Tables?|Notes?|panels?|Reviewer|[Cc]oncerns?|Comments?|[Ii]tems?|[Pp]aragraphs?|Mechanism|Step|Phase|[Ll]imitations?|[Ss]ections?|[Aa]dditional file|[Pp]ages?|[Ee]quation|Eq\.|[Pp]oints?|Supplementary Notes?)"
THRESH = {0.05, 0.1, 0.25, 0.3, 0.5, 0.75, 0.9, 0.95, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 0.2, 0.8, 0.6, 0.4}
def struct_reason(tok, num, sec, before, after, sent, fname, is_fig):
    v = None
    try: v = float(num.replace(",", ""))
    except Exception: pass
    if sec in ("References", "Title page"): return f"{sec.lower()}"
    b0 = before[-60:]; a0 = after[:40]
    if re.match(r"\s*threads?\b", a0): return "software parameter (threads)"
    if re.search(r"Build\s*$", b0): return "assembly build number"
    if re.match(r"-dpi\b", a0): return "figure format"
    if re.search(r"[–—\-]{1,2}[A-Za-z]\w*\s+$", b0): return "command-line parameter"
    if sec == "Methods" and re.search(r"quality", sent) and re.search(r"were\s*$|of\s*$|≥\s*$", b0): return "Methods parameter (quality filter)"
    if re.search(r"\btop\s*$", b0) and re.search(r"display|Heatmap|heatmap", sent): return "display choice"
    if re.search(r"p\s*<\s*$", b0) and v in (0.05, 0.01, 0.001) and "*" in sent: return "significance legend"
    if re.search(r"(?:Δ|line at|from)\s*=?\s*$", b0) and v in (0.0, 0.5) and re.search(r"dashed|line|deviation from", sent): return "reference line / null value"
    if re.match(r"^\s*\d+\.\s*$", sent) or re.match(r"^\s*" + re.escape(num) + r"\.\s{1,3}\S", sent): return "list ordinal"
    if re.search(r"\b\d{4};\s*\d+", sent) and re.search(r"\d{4};\s*[\d()]*$|\(\s*$|\d{4}\s*$|;\s*$", b0): return "bibliographic citation"
    if re.search(r"[Rr]eferences?\s*$", b0): return "reference number"
    if re.match(r"%\s*confidence", tok[-1:] + a0) or (tok.endswith("%") and re.match(r"\s*confidence", a0)): return "confidence level"
    b = before[-60:]; a = after[:40]
    if re.search(r"\[[\d,–\- ]*$", before) and re.match(r"[\d,–\- ]*\]", after): return "citation number"
    if re.search(LABEL + r"\s*(?:S?\d+[A-Za-z]?\s*(?:[–,]|and|or|to)\s*)*S?$", b): return "figure/table/note/panel or letter label"
    if re.search(r"(?:chromosomes?|chr\.?|Chr)\s*(?:\d+\s*(?:[–,/]|and|or)\s*)*$", b): return "chromosome number"
    if re.search(r"\bseed\s*$", b): return "random seed (Methods)"
    if re.search(r"(?:Ensembl|[Rr]elease|REL-?)\s*$", b): return "release number"
    if re.search(SOFTWARE + r"\s*\(?\s*(?:version\s*|v)?$", b): return "software version"
    if re.search(r"(?:(?:^|[\s(])--?[A-Za-z][\w-]*|\b[a-z]+[A-Z]\w*)\s*[= ]\s*$", b): return "command-line parameter"
    if re.search(r"Verdict rules|verdict rules|was called|were called|is called", sent) and sec == "Methods": return "identity verdict rule (fixed in Methods)"
    if sec == "Methods" and re.search(r"(?:≥|≤|>|<|at least|at most|between)\s*$", b): return "threshold fixed in Methods"
    if re.search(r"\bage\s*$", b) and re.search(r"GEO|recorded", sent): return "GEO sample attribute"
    if re.search(r"(?:≥|>)\s*$", b) and v in (80.0, 85.0, 90.0, 0.8, 0.85) and re.search(r"residual|homozygous|sites", sent): return "definition threshold (Methods/Note)"
    if re.match(r"\s*\)\s*$", a) and re.search(r"\(\s*$", b) and v is not None and v < 10 and not tok.count("."): pass
    if v is not None and "." not in num and "," not in num and 1950 <= v <= 2030: return "year"
    if re.match(r"^\s*\d+\.\s{1,3}\S", sent) and sent.lstrip().startswith(num + "."): return "list ordinal"
    if re.match(r"\s*(?:-|‑|\s)?bp\b|\s*-bp\b", a) and (re.search(r"read|Read|single-end|paired|Paired|library", sent) or sec == "Table 1"): return "read length (design descriptor)"
    if re.match(r"\s*(?:bp|nt)\b", a) and sec == "Table 1": return "read length (design descriptor)"
    if v is not None and v in THRESH and re.search(r"(?:≥|≤|>|<|±|\|Δ\|\s*[>≤<≥]|at least|threshold of|cut-?off of|FDR\s*[<≤])\s*$", b): return "threshold fixed in Methods"
    if v is not None and v in THRESH and re.match(r"\s*log₂", a) and re.search(r"(?:threshold|≥|bias of|exceed|criterion|cut)", sent): return "threshold fixed in Methods"
    if v == 0.5 and re.search(r"expect|parity|equal|balanced", sent): return "design constant (expected F1 allele fraction)"
    if v in (0.5, 0.50) and "%" not in tok and re.search(r"0\.50?\s*log₂|bias\s*[≥>]", sent) : return "threshold fixed in Methods"
    if v == 255 and re.search(r"MAPQ|mapping quality|uniquely", sent): return "STAR MAPQ (parameter)"
    if v == 20 and re.search(r"depth|reads per site|≥ 20|total count", sent): return "site depth threshold (Methods)"
    if v == 10 and "%" in tok and re.search(r"top|bottom|decile|subset", sent): return "decile/subset definition (Methods)"
    if v is not None and re.search(r"(?:CO₂|°C|μg|µg|mg|ml|μl|nM|mM)", a[:4]): return "unit constant"
    if re.search(r"(?:Tel|Fax|phone)[:.]?\s*\+?[\d\- ]*$", b) or re.search(r"Aomori|Towada|Higashi", sent): return "address/phone"
    return None

def scored(keys, sent):
    words = set(re.findall(r"[a-z0-9]+", sent.lower()))
    return sorted(keys, key=lambda k: (-len(set(k.lower().split("_")) & words), len(k), k))[:3]

def classify(tok, kind, m, sec, before, after, sent, fname, is_fig):
    if kind == "plain":
        num = m.group("num"); pct = m.group("pct")
        r = struct_reason(tok, num, sec, before, after, sent, fname, is_fig)
        if r: return "STRUCT", r
        d = len(num.split(".")[1]) if "." in num else 0
        s = num.replace(",", "")
        cand = {s}
        ks = KEY_IDX.get((d, s), set())
        if ks: return "KEY", ";".join(scored(ks, sent))
        xs = XL_IDX.get((d, s), set())
        if xs: return "XLSX", ";".join(sorted(xs, key=len)[:3])
        return "UNMAPPED", ""
    if kind in ("sci", "enot"):
        if kind == "sci": mant = m.group("sm").replace("−", "-").lstrip("+-"); e = int(m.group("se").translate(SUP))
        else: mant, e = tok.split("e"); e = int(e)
        d = len(mant.split(".")[1]) if "." in mant else 0
        mv = float(mant);
        # normalise e.g. 0.0e+00
        if mv == 0:
            ks = KEY_ZERO_P
            return ("KEY", ";".join(scored(ks, sent)) + " (printed 0.0e+00: exact value underflows)") if ks else ("UNMAPPED", "")
        ks = KEY_SCI.get((d, f"{mv:.{d}f}", e), set())
        if ks: return "KEY", ";".join(scored(ks, sent))
        return "UNMAPPED", ""
    if kind == "pow":
        e = int(m.group("pe").translate(SUP))
        if is_fig: return "STRUCT", "axis tick (log scale)"
        if re.search(r"[<≤]\s*$", before):
            ks = [k for k, x in KEYS if x < 10 ** e and ("_p" in k or k.endswith("_p") or "pval" in k)]
            if ks: return "KEY", ";".join(scored(ks, sent)) + f" (bound < 10^{e})"
        return "UNMAPPED", ""
    return "UNMAPPED", ""

# ---------------------------------------------------------------- documents ----------------------------------------------------------
MAIN_H = {"Abstract": "Abstract", "Keywords": "Keywords", "Background": "Background", "Methods": "Methods", "Results": "Results",
          "Discussion": "Discussion", "Conclusions": "Conclusions", "List of abbreviations": "List of abbreviations",
          "Declarations": "Declarations", "References": "References", "Figure legends": "Figure legends"}
def sections_main(lines):
    sec = "Title page"; out = []
    for ln in lines:
        t = ln.strip()
        if sec == "Abstract" and t in ("Background", "Methods", "Results", "Conclusions"): out.append((sec, None)); continue
        if t in MAIN_H and not (sec == "Abstract" and t != "Keywords"): sec = MAIN_H[t]; out.append((sec, None)); continue
        if t.startswith("Table 1."): sec = "Table 1"
        out.append((sec, ln))
    return out
def sections_generic(lines, first):
    sec = first; out = []
    for ln in lines:
        t = ln.strip()
        if t and len(t) < 90 and not re.search(r"[.:;,)]$", t) and not t.startswith(("-", "|", "+")) and not re.match(r"^[\d\s.,−\-+%]+$", t):
            sec = t[:80]; out.append((sec, None)); continue
        out.append((sec, ln))
    return out
def split_sent(par):
    spans = []; start = 0
    for m in re.finditer(r"(?<!Fig)(?<!Figs)(?<!al)(?<!vs)(?<!e\.g)(?<!i\.e)(?<!ca)(?<!no)[.;](?=\s+[A-Z(*])", par):
        spans.append((start, m.end())); start = m.end()
    spans.append((start, len(par)))
    return [(a, b) for a, b in spans if par[a:b].strip()]

rows = []
TEXT_NUMS = collections.defaultdict(set)    # normalised numeric string -> text files
def norm_num(s): return s.replace("−", "-").replace(",", "").lstrip("+-±")
def audit_text(fname, pairs, is_fig=False):
    for sec, ln in pairs:
        if ln is None or not ln.strip(): continue
        clean = STRIP.sub(blank, ln)
        for a, b in (split_sent(clean) if not is_fig else [(0, len(clean))]):
            sent = clean[a:b]
            for m in TOK.finditer(sent):
                kind = m.lastgroup if m.lastgroup in ("sci", "enot", "pow", "plain") else next(k for k in ("sci", "enot", "pow", "plain") if m.group(k))
                tok = m.group(0).strip()
                before, after = sent[:m.start()], sent[m.end():]
                if re.search(r"[A-Za-z]/$", before): continue   # strain identifiers such as C57BL/6
                if is_fig:
                    cls, ref = fig_class(tok, kind, m, sec, before, after, sent, fname)
                else:
                    cls, ref = classify(tok, kind, m, sec, before, after, sent, fname, False)
                    if kind == "plain": TEXT_NUMS[norm_num(m.group("num")) + m.group("pct")].add(fname)
                    else: TEXT_NUMS[norm_num(tok)].add(fname)
                note = ""
                if tok in ("0.333", "0.386") and re.search(r"0\.333|submitted|not 0\.333|naive B6 mapping", sent):
                    cls, ref, note = "UNMAPPED", "", "submitted-version (v1) F1 value quoted"
                elif cls == "KEY":
                    words = set(re.findall(r"[a-z0-9]+", sent.lower()))
                    best = ref.split(";")[0].split(" ")[0]
                    if not (set(best.lower().split("_")) & words): note = "weak: no word of the key name in the sentence"
                if cls != "STRUCT" and (tok.startswith("±") or re.search(r"(?:within|fewer than|less than|more than|at most|at least|below|above|up to|<|>|≤|≥)\s*[±]?$", before[-16:])):
                    note = (note + "; " if note else "") + "bound: printed value is a bound — PASS only if the exact value satisfies it (M-2)"
                rows.append([fname, sec, re.sub(r"\s+", " ", sent.strip())[:300], tok, cls, ref, note])

TICK_VALUES = {}
STEPS = (0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 10000)
def tick_sets(txt):
    """values of standalone numbers (separated by >= 2 spaces / line ends) in a figure, keyed by decimal count"""
    vals = collections.defaultdict(set)
    for ln in txt.splitlines():
        for m in re.finditer(r"(?:(?<=\s\s)|^)([−\-]?\d+(?:\.\d+)?)(?=\s\s|$)", ln):
            t = m.group(1); d = len(t.split(".")[1]) if "." in t else 0
            vals[d].add(round(float(t.replace("−", "-")), 6))
    return vals
def is_tick(v, d, vals):
    S = vals.get(d, set())
    for st in STEPS:
        a, b2, c, e = round(v - st, 6), round(v + st, 6), round(v - 2 * st, 6), round(v + 2 * st, 6)
        if (a in S and b2 in S) or (b2 in S and e in S) or (a in S and c in S): return True
    return False
def fig_class(tok, kind, m, sec, before, after, sent, fname):
    if kind == "pow": return "STRUCT", "axis tick (log scale)"
    if kind == "plain" and re.search(r"clipped to \[[^\]]*$", before): return "STRUCT", "display range (clipping)"
    if kind == "plain" and re.search(r"\b(?:heart|rep|replicate|animal)\s*$", before) and re.match(r"\s*\(GSM", after): return "STRUCT", "sample label"
    if kind == "plain" and re.match(r"-Mb\b", after): return "STRUCT", "bin size (display)"
    if kind == "plain" and re.search(r"\btop\s*$", before): return "STRUCT", "display choice"
    if kind == "plain" and m.group("sign") != "+" and not re.search(r"(?:=|mean|ρ|p|n|of|\()\s*$", before[-8:]):
        num = m.group("num"); d = len(num.split(".")[1]) if "." in num else 0
        v = float(num.replace(",", "")) * (-1 if m.group("sign") in ("−", "-") else 1)
        if "," not in num and is_tick(round(v, 6), d, FIGVALS) and re.search(r"(^|\s{2,})$", before) and re.match(r"(\s{2,}|$)", after):
            return "STRUCT", "axis/colour-bar tick (arithmetic progression)"
    standalone = bool(re.search(r"(^|\s{2,})$", before)) and bool(re.match(r"(\s{2,}|$)", after))
    lab = re.search(r"(?:=|mean|n\s*=|of|ρ|p\s*[=<]|SE|±|expected|median|r\s*=)\s*$", before[-14:]) or re.match(r"\s*(?:of|genes|sites|catalog|SNPs|animals|samples|%)", after)
    if kind == "plain" and not lab and m.group("sign") in ("", "−", "-"):
        num = m.group("num"); v = float(num.replace(",", ""))
        d = len(num.split(".")[1]) if "." in num else 0
        # standalone round values (ticks of linear axes and colour bars)
        if standalone and (d <= 1 or (d == 2 and round(v * 100) % 5 == 0)):
            return "STRUCT", "axis/colour-bar tick"
        # a standalone run of round numbers on one line (x-axis tick rows)
        nums = re.findall(r"(?<![\w.])[−\-]?\d+(?:\.\d+)?(?![\w.])", sent)
        if len(nums) >= 3 and re.fullmatch(r"[\s\d.−\-]+", sent.strip()) and d <= 2:
            return "STRUCT", "axis tick row"
    if kind == "plain" and re.search(r"(?:\bpanel|^\s*)$", before) and re.match(r"\s*$", after) and len(tok) == 1: return "STRUCT", "panel"
    if kind == "plain" and re.search(r"10\s*$", before[-4:]): pass
    r = struct_reason(tok, m.group("num") if kind == "plain" else tok, sec, before, after, sent, fname, True) if kind == "plain" else None
    if r: return "STRUCT", r
    return classify(tok, kind, m, sec, before, after, sent, fname, True)

for fn in ("LAR_main_text_v2.plain.txt",):
    audit_text(fn, sections_main(open(f"{SUB}/{fn}", encoding="utf-8").read().splitlines()))
for fn, first in (("LAR_additional_file_1.plain.txt", "Additional file 1"), ("response_to_reviewers_v2.plain.txt", "Response to Reviewers"), ("cover_letter_v2.plain.txt", "Cover letter")):
    audit_text(fn, sections_generic(open(f"{SUB}/{fn}", encoding="utf-8").read().splitlines(), first))
for pdf in sorted(glob.glob(f"{SUB}/figures/*.pdf")):
    txt = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, text=True).stdout
    name = os.path.basename(pdf)
    FIGVALS = tick_sets(txt)
    audit_text(f"figures/{name}", [(name[:-4], ln) for ln in txt.splitlines()], is_fig=True)

# figure-vs-text note
for r in rows:
    if r[0].startswith("figures/") and r[4] != "STRUCT":
        t = r[3]
        t2 = t.replace("e-", "e-")
        k = norm_num(re.sub(r"%$", "", t)) + ("%" if t.endswith("%") else "")
        hits = sorted(TEXT_NUMS.get(k, set()) | TEXT_NUMS.get(norm_num(t), set()))
        if "e" in t and re.fullmatch(r"\d+(\.\d+)?e[-+]\d+", t):
            mnt, e = t.split("e"); hits = sorted(TEXT_NUMS.get(f"{mnt} × 10{str(int(e)).translate(str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻'))}", set()))
            if not hits:
                for kk, fs in TEXT_NUMS.items():
                    if kk.startswith(mnt) and "×" in kk and kk.endswith(str(int(e)).translate(str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻'))): hits = sorted(fs)
        r[6] = (r[6] + "; " if r[6] else "") + ("text: " + ",".join(h.split(".")[0] for h in hits)) if hits else "no text counterpart"

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8", newline="") as o:
    w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["file", "section", "sentence", "token", "class", "ref", "note"]); w.writerows(rows)
c = collections.Counter(r[4] for r in rows); print("tokens", len(rows), dict(c))
per = collections.defaultdict(collections.Counter)
for r in rows: per[r[0]][r[4]] += 1
for f, cc in per.items(): print(f"  {f:55s} {dict(cc)}")
