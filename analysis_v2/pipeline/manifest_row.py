#!/usr/bin/env python3
"""Print one manifest row as shell variable assignments. Usage: manifest_row.py <sample_id>"""
import csv, sys, shlex
PROJ="/usr/local/jupyter/FVB_bias"
sid=sys.argv[1]
for r in csv.DictReader(open(f"{PROJ}/revision_v2/manifest/manifest.tsv"), delimiter="\t"):
    if r["sample_id"]==sid:
        for k,v in r.items(): print(f"M_{k}={shlex.quote(v)}")
        sys.exit(0)
print(f"echo 'sample_id {sid} not in manifest' >&2; exit 1")
