#!/usr/bin/env python3
"""Phase H-4: convert [CITE Author Year, ...] placeholders to numbered citations [n] in order of first appearance across the
main text (Background → Methods → Results → Discussion → Conclusions → figure legends → Table 1 → Additional file), and write
references.md renumbered in Vancouver style with DOI. Software/database entries are given numbers like any other reference.
Usage: number_citations.py <section files in citation order...> --refs <references_source.md> --out-refs <references.md>
The source reference file is manuscript/revised_v2/references_source.md: one entry per line, "KEY | formatted reference", where
KEY is the [CITE] key (e.g. "Keane 2011") or a software key used in the text (e.g. "STAR"). Entries never cited are reported."""
import re, sys
args=sys.argv[1:]; refs_src=args[args.index("--refs")+1]; out_refs=args[args.index("--out-refs")+1]; files=[a for a in args if not a.startswith("--") and a not in (refs_src,out_refs)]
entries={}
for l in open(refs_src, encoding="utf-8"):
    if "|" not in l or l.startswith("#"): continue
    k,v=l.split("|",1); entries[k.strip()]=v.strip()
order=[]; missing=set()
def repl(m):
    keys=[k.strip() for k in m.group(1).split(",")]; nums=[]
    for k in keys:
        if k not in entries: missing.add(k); return m.group(0)
        if k not in order: order.append(k)
        nums.append(order.index(k)+1)
    nums=sorted(set(nums)); # compress runs: 1,2,3 -> 1–3
    parts=[]; i=0
    while i<len(nums):
        j=i
        while j+1<len(nums) and nums[j+1]==nums[j]+1: j+=1
        parts.append(f"{nums[i]}–{nums[j]}" if j>i+1 else (f"{nums[i]},{nums[j]}" if j==i+1 else f"{nums[i]}")); i=j+1
    return "["+",".join(parts)+"]"
for f in files:
    s=open(f,encoding="utf-8").read(); s2=re.sub(r"\[CITE ([^\]]+)\]", repl, s)
    if s2!=s: open(f,"w",encoding="utf-8").write(s2)
with open(out_refs,"w",encoding="utf-8") as o:
    o.write("# References\n\n")
    for i,k in enumerate(order,1): o.write(f"{i}. {entries[k]}\n\n")
print(f"{len(order)} references numbered; cited keys in order: {order}")
uncited=[k for k in entries if k not in order]; print("uncited entries:", uncited); print("keys without entry:", sorted(missing))
sys.exit(1 if missing else 0)
