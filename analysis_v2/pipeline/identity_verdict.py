#!/usr/bin/env python3
"""Strain-identity verdict for one sample (same metrics/rules as revision_v2/scripts/strain_identity_qc.py, fixed before Phase A).
Usage: identity_verdict.py <ase.tsv> <manifest_strain>  -> prints a TSV line; exit 1 if the verdict does not match the expected class."""
import sys
EXPECT={"FVB/NJ":"FVB/NJ pure","C57BL/6J":"B6 pure","B6.FVB-Tg(MMTV-PyVT)634Mul/LellJ":"B6 pure"}
def verdict(m,hi,mid,lo):
    if m>0.90 and hi>0.85: return "FVB/NJ pure"
    if m<0.10 and lo>0.85: return "B6 pure"
    if mid>0.50 and 0.25<=m<=0.55: return "F1 (B6 x FVB)"
    if hi>0.15 and mid>0.15 and lo>0.15: return "segregating mixed background"
    if hi>0.25 and lo>0.25 and mid<0.15: return "other inbred strain"
    return "needs review"
afs=[]
with open(sys.argv[1]) as f:
    h=f.readline().rstrip("\n").split("\t"); i={k:j for j,k in enumerate(h)}
    for l in f:
        r=l.rstrip("\n").split("\t"); ref,alt=int(r[i["refCount"]]),int(r[i["altCount"]]); t=ref+alt
        if t>=20: afs.append(alt/t)
n=len(afs); m=sum(afs)/n if n else float("nan"); hi=sum(a>0.9 for a in afs)/n if n else 0; lo=sum(a<0.1 for a in afs)/n if n else 0; mid=1-hi-lo if n else 0
v=verdict(m,hi,mid,lo) if n else "needs review (no sites)"
exp=EXPECT.get(sys.argv[2],"?")
print("n_sites\tmean_af\tfrac_hi\tfrac_mid\tfrac_lo\tqc_verdict\texpected\tmatch")
ok= v==exp
print(f"{n}\t{m:.3f}\t{hi:.3f}\t{mid:.3f}\t{lo:.3f}\t{v}\t{exp}\t{'Y' if ok else 'N'}")
sys.exit(0 if ok else 1)
