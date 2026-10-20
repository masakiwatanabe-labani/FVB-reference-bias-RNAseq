#!/usr/bin/env python3
"""Strain-identity QC from ASEReadCounter output at FVB/NJ homozygous-alt exonic SNP sites.
Metrics per sample: n_sites, mean_af (alt fraction), frac_hi (AF>0.9), frac_mid (0.1-0.9), frac_lo (AF<0.1).
Verdict thresholds are fixed in advance (revision instructions v2, A-2) and are NOT tuned to the data.
Usage: strain_identity_qc.py <ase_dir> <out.tsv> [min_depth]
"""
import sys, os, glob

def verdict(mean_af, hi, mid, lo):
    if mean_af > 0.90 and hi > 0.85: return "FVB/NJ pure"
    if mean_af < 0.10 and lo > 0.85: return "B6 pure"
    if mid > 0.50 and 0.25 <= mean_af <= 0.55: return "F1 (B6 x FVB)"
    if hi > 0.15 and mid > 0.15 and lo > 0.15: return "segregating mixed background"
    if hi > 0.25 and lo > 0.25 and mid < 0.15: return "other inbred strain"
    return "needs review"

def summarise(path, min_depth=20):
    afs = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        i = {h: k for k, h in enumerate(hdr)}
        for line in f:
            r = line.rstrip("\n").split("\t")
            ref, alt = int(r[i["refCount"]]), int(r[i["altCount"]])
            tot = ref + alt
            if tot < min_depth: continue
            afs.append(alt / tot)
    n = len(afs)
    if n == 0:
        return n, float("nan"), float("nan"), float("nan"), float("nan"), "needs review (no sites)"
    mean_af = sum(afs) / n
    hi = sum(a > 0.9 for a in afs) / n
    lo = sum(a < 0.1 for a in afs) / n
    mid = 1 - hi - lo
    return n, mean_af, hi, mid, lo, verdict(mean_af, hi, mid, lo)

if __name__ == "__main__":
    ase_dir, out = sys.argv[1], sys.argv[2]
    md = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    with open(out, "w") as o:
        o.write("SRR\tn_sites\tmean_af\tfrac_hi\tfrac_mid\tfrac_lo\tqc_verdict\n")
        for p in sorted(glob.glob(os.path.join(ase_dir, "*.tsv"))):
            srr = os.path.basename(p)[:-4]
            n, m, hi, mid, lo, v = summarise(p, md)
            o.write(f"{srr}\t{n}\t{m:.3f}\t{hi:.3f}\t{mid:.3f}\t{lo:.3f}\t{v}\n")
    print(open(out).read())
