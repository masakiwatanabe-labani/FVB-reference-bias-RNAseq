#!/usr/bin/env python3
"""Phase M-3 D3: double-rounding scan. For every numbers.tsv key printed in the text («key|fmt» in build_phase_j.py, the donor files and
the section files; f-string keys such as rps2mm_{g}_… expanded over the four Rps2 genes), compare the printed value obtained from the stored
value with the value obtained by rounding the full-precision value once (generated/numbers_fullprecision.tsv). Any difference is a
double-rounding defect. Usage: double_rounding_scan.py [numbers.tsv]   (default: the current numbers.tsv). Exit 1 if any defect."""
import pandas as pd, re, glob, sys
P = "/usr/local/jupyter/FVB_bias"; M = f"{P}/manuscript/revised_v2"
fp = pd.read_csv(f"{M}/generated/numbers_fullprecision.tsv", sep="\t", dtype=str, keep_default_na=False).drop_duplicates("key", keep="last").set_index("key").value
used = {}
for f in [f"{P}/revision_v2/scripts/build_phase_j.py"] + glob.glob(f"{P}/revision_v2/audit/pre_numbering/*.md") + glob.glob(f"{M}/*.md"):
    for k, fm in re.findall(r"«([^|»]+)\|([^»]+)»", open(f, encoding="utf-8").read()):
        for kk in ([k.replace("{g}", g) for g in ("Rps2", "Rps2_ps10", "Rps2_ps13", "Rps2_ps6")] if "{g}" in k else [k]): used.setdefault(kk, set()).add(fm)
st = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else f"{M}/generated/numbers.tsv", sep="\t", dtype=str, keep_default_na=False).drop_duplicates("key", keep="last").set_index("key").value
out = []
for k in sorted(used):
    if k not in st.index or k not in fp.index: continue
    try: x6, xf = float(st[k]), float(fp[k])
    except ValueError: continue
    for fm in used[k]:
        m = re.fullmatch(r"[+]?\.(\d)f", fm) or re.fullmatch(r"pct(\d)", fm)
        if m: d = int(m.group(1)); a, b = f"{x6:.{d}f}", f"{xf:.{d}f}"
        elif fm in ("d", ",d"): a, b = str(int(round(x6))), str(int(round(xf)))
        elif fm in ("abs2", "abs3"): d = int(fm[-1]); a, b = f"{abs(x6):.{d}f}", f"{abs(xf):.{d}f}"
        else: continue
        if a != b: out.append((k, fm, st[k], fp[k], a, b))
print(f"keys printed in the text: {len(used)}; checked against full precision: {sum(1 for k in used if k in st.index)}; double-rounding defects: {len(out)}")
for r in out: print("  ", r)
sys.exit(1 if out else 0)
