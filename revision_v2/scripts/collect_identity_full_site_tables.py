#!/usr/bin/env python3
"""Phase M-3 D1b: per-site allele-count tables for Supplementary Fig. S1, collected into analysis_v2/results/identity/full_site_tables/.
Reference animals by rule (no hand choice): among the animals evaluated by the v2 pipeline at the full 203,510-site set (a
<sample>.ase.tsv containing sites outside the 20,351-site subset), whose identity verdict matches the recorded strain, take for each strain the
animal with the most sites of depth >= 20 (ties: first sample_id in sort order). The four GSE55788 heart animals are the subject of the
figure (verdict "needs review") and are not candidates. Writes candidates.tsv (every candidate, its site count, the choice) and copies
the six per-site tables. Usage: collect_identity_full_site_tables.py"""
import glob, gzip, os, shutil, pandas as pd
P = "/usr/local/jupyter/FVB_bias"; D = f"{P}/analysis_v2/data"; OUT = f"{P}/analysis_v2/results/identity/full_site_tables"; os.makedirs(OUT, exist_ok=True)
man = pd.read_csv(f"{P}/revision_v2/manifest/manifest.tsv", sep="\t", dtype=str).set_index("sample_id")
sub = set()
for l in gzip.open(f"{P}/analysis_v2/qc/fvbnj_homalt_exonic.hetGT.sub10pct.seed20261016.vcf.gz", "rt"):
    if l[0] != "#": c_, p_ = l.split("\t")[:2]; sub.add((c_, int(p_)))
EXPECT = {"C57BL/6J": "B6 pure", "FVB/NJ": "FVB/NJ pure"}
rows = []
for f in sorted(glob.glob(f"{D}/*/*/qc/*.ase.tsv")):
    s = os.path.basename(f)[:-8]
    if s.startswith("heart_") or s not in man.index: continue
    a = pd.read_csv(f, sep="\t", usecols=["contig", "position", "totalCount"], dtype={"contig": str}); idv = pd.read_csv(os.path.join(os.path.dirname(f), "identity.tsv"), sep="\t").iloc[0]
    if all((c_, p_) in sub for c_, p_ in zip(a.contig, a.position)): continue   # evaluated at the subset only
    st = man.strain[s]
    if st not in EXPECT or idv.qc_verdict != EXPECT[st]: continue
    rows.append(dict(sample_id=s, strain=st, tissue=man.tissue[s], GSE=man.GSE[s], GSM=man.GSM[s], n_sites_depth_ge_20=int((a.totalCount >= 20).sum()), identity_verdict=idv.qc_verdict, source=os.path.relpath(f, P)))
c = pd.DataFrame(rows).sort_values(["strain", "n_sites_depth_ge_20", "sample_id"], ascending=[True, False, True])
c["chosen"] = ""
for st in EXPECT: c.loc[c[c.strain == st].index[0], "chosen"] = "Y"
c.to_csv(f"{OUT}/candidates.tsv", sep="\t", index=False); print(c.to_string(index=False))
copy = [f"{D}/heart/heart_B6_C57_M_12wk_mRNA_{i}/qc/heart_B6_C57_M_12wk_mRNA_{i}.ase.tsv" for i in range(1, 5)] + [f"{P}/{x}" for x in c[c.chosen == "Y"].source]
for f in copy:
    shutil.copy2(f, OUT); print(f"{os.path.basename(f)}\t{os.path.getsize(f)/1e6:.2f} MB")
