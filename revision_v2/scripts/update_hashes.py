#!/usr/bin/env python3
"""Phase W (2026-10-26): refresh the hash lists after a regeneration, instead of editing them by hand.
For every listed path that matches one of the PATTERNS (fnmatch on the path as written in the list), the hash is recomputed; paths in
--add are appended when absent. Lists and their base directories:
  analysis_v2/results/RESULTS.md5 (md5, paths relative to analysis_v2/results, written "./x"), analysis_v2/figures/FIGURES.md5 (md5,
  relative to analysis_v2/figures), revision_v2/freeze/inputs.sha256 and deliverables.sha256 (sha256, relative to the project).
Prints every changed or added line. Usage: update_hashes.py PATTERN [PATTERN ...] [--add LIST:PATH ...]"""
import sys, os, hashlib, fnmatch
P = "/usr/local/jupyter/FVB_bias"
LISTS = {"RESULTS": (f"{P}/analysis_v2/results/RESULTS.md5", f"{P}/analysis_v2/results", "md5"),
         "FIGURES": (f"{P}/analysis_v2/figures/FIGURES.md5", f"{P}/analysis_v2/figures", "md5"),
         "inputs": (f"{P}/revision_v2/freeze/inputs.sha256", P, "sha256"),
         "deliverables": (f"{P}/revision_v2/freeze/deliverables.sha256", P, "sha256")}
args = sys.argv[1:]; adds = []
if "--add" in args: i = args.index("--add"); adds = args[i + 1:]; args = args[:i]
def h(path, algo):
    x = hashlib.new(algo)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): x.update(b)
    return x.hexdigest()
for name, (lst, base, algo) in LISTS.items():
    lines = open(lst).read().splitlines(); out = []; changed = 0
    for l in lines:
        hv, p = l.split(None, 1)
        if any(fnmatch.fnmatch(p, pat) or fnmatch.fnmatch(p.lstrip("./"), pat) for pat in args):
            nv = h(os.path.join(base, p), algo)
            if nv != hv: print(f"{name}: {p}  {hv[:8]} -> {nv[:8]}"); changed += 1
            l = f"{nv}  {p}"
        out.append(l)
    present = {l.split(None, 1)[1] for l in out}
    for a in adds:
        ln, p = a.split(":", 1)
        if ln == name and p not in present:
            out.append(f"{h(os.path.join(base, p), algo)}  {p}"); print(f"{name}: added {p}"); changed += 1
    if changed: open(lst, "w").write("\n".join(out) + "\n")
