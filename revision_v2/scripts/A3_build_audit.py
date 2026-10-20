#!/usr/bin/env python3
"""Phase A-3: build revision_v2/audit/sample_audit.tsv from
  - old_runs_metadata.tsv (ENA + BioSample + GEO SOFT, per SRR)
  - analysis_v2/qc/strain_identity_qc.tsv (ASEReadCounter-based verdicts)
  - the manuscript's stated design (hard-coded here from submission/manuscript_combined.md L35-37, Table 1)
"""
import collections
PROJ = "/usr/local/jupyter/FVB_bias"; AUD = f"{PROJ}/revision_v2/audit"; QC = f"{PROJ}/analysis_v2/qc"

# Manuscript claims (strain, tissue, expected QC verdict) per GSE and old arm
CLAIM = {
 ("GSE123875","FVB"): ("FVB/NJ", "mesenteric adipose", "FVB/NJ pure"),
 ("GSE123875","B6"):  ("C57BL/6J", "mesenteric adipose", "B6 pure"),
 ("GSE123893","FVB"): ("FVB/NJ", "liver", "FVB/NJ pure"),
 ("GSE123893","B6"):  ("C57BL/6J", "liver", "B6 pure"),
 ("GSE123894","FVB"): ("FVB/NJ", "kidney", "FVB/NJ pure"),
 ("GSE123894","B6"):  ("C57BL/6J", "kidney", "B6 pure"),
 ("GSE135230","FVB"): ("FVB/N", "cardiac muscle", "FVB/NJ pure"),
 ("GSE135230","B6"):  ("C57BL/6J", "cardiac muscle", "B6 pure"),
 ("GSE175625","FVB"): ("FVB/N", "mammary tumour (PyVT)", "FVB/NJ pure"),
 ("GSE175625","B6"):  ("C57BL/6J", "mammary tumour (PyVT)", "B6 pure"),
 ("GSE200632","F1"):  ("C57BL/6 x FVB/N F1", "primary keratinocytes", "F1 (B6 x FVB)"),
 ("GSE293440","F1"):  ("C57BL/6J x FVB/NJ F1", "frontal cortex", "F1 (B6 x FVB)"),
}
TISSUE_SYN = {"mesenteric adipose": ["adipose"], "liver": ["liver"], "kidney": ["kidney"],
              "cardiac muscle": ["heart", "cardiac"], "mammary tumour (PyVT)": ["mammary", "tumor", "tumour"],
              "primary keratinocytes": ["keratinocyte"], "frontal cortex": ["frontal", "cortex"]}

meta = {}
hdr = None
for l in open(f"{AUD}/old_runs_metadata.tsv").read().splitlines():
    r = l.split("\t")
    if hdr is None: hdr = r; continue
    meta[r[0]] = dict(zip(hdr, r))
qc = {}
for l in open(f"{QC}/strain_identity_qc.tsv").read().splitlines()[1:]:
    r = l.split("\t"); qc[r[0]] = r

cols = ["GSE","GSM","SRR","BioSample","geo_strain","geo_genotype","geo_tissue","geo_treatment","geo_sex","individual_id",
        "n_technical_runs","old_group","n_sites","mean_af","frac_hi","frac_mid","frac_lo","qc_verdict","matches_manuscript","evidence_files"]
out = open(f"{AUD}/sample_audit.tsv", "w"); out.write("\t".join(cols) + "\n")
for srr, m in meta.items():
    gse, grp = m["GSE"], m["old_group"]
    cs, ct, cv = CLAIM[(gse, grp)]
    q = qc.get(srr, [srr, "NA", "NA", "NA", "NA", "NA", "QC not run"])
    strain = m["geo_strain"] or m["biosample_strain"]
    tissue = m["geo_tissue"] or m["geo_source_name"]
    geno = m["geo_genotype"]
    reasons = []
    # strain: manuscript strain string must be the stated strain, not merely a substring of a background name
    s_norm = strain.replace(" ", "").lower(); cs_norm = cs.replace(" ", "").lower()
    if grp == "F1":
        strain_ok = ("fvb" in s_norm and ("c57" in s_norm or "b6" in s_norm)) or "f1" in s_norm
    else:
        strain_ok = s_norm == cs_norm or (cs_norm == "fvb/n" and s_norm in ("fvb/n", "fvb/nj")) or (cs_norm == "c57bl/6j" and s_norm in ("c57bl/6j", "c57bl/6"))
    if not strain_ok: reasons.append(f"strain: manuscript '{cs}' vs GEO '{strain}'")
    tissue_ok = any(k in tissue.lower() for k in TISSUE_SYN[ct])
    if not tissue_ok: reasons.append(f"tissue: manuscript '{ct}' vs GEO '{tissue}'")
    if q[6] != cv: reasons.append(f"QC: expected '{cv}' got '{q[6]}'")
    if int(m["n_technical_runs"]) > 1: reasons.append(f"technical run ({m['n_technical_runs']} runs share BioSample)")
    match = "Y" if not reasons else "N: " + "; ".join(reasons)
    ev = f"revision_v2/audit/raw/biosample_{m['BioSample']}_*.xml; revision_v2/audit/raw/{gse}_family.soft_*.gz; revision_v2/audit/raw/ena_filereport_{srr}_*.tsv; analysis_v2/qc/ase/{srr}.tsv"
    row = [gse, m["GSM"], srr, m["BioSample"], strain, geno, tissue, m["geo_treatment"], m["geo_sex"], m["geo_title"] or m["sample_title_ena"],
           m["n_technical_runs"], grp, q[1], q[2], q[3], q[4], q[5], q[6], match, ev]
    out.write("\t".join(row) + "\n")
out.close()
print(open(f"{AUD}/sample_audit.tsv").read())
