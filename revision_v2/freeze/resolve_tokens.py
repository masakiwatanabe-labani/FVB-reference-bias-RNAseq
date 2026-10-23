#!/usr/bin/env python3
"""Phase M-1/M-2 link (rerun after token_audit.py): annotate token_audit.tsv with the M-2 claim that resolves each UNMAPPED token (and the coincidental XLSX
match of Fig. 2C), and the claim status from recompute_report.tsv. Adds columns resolved_by, resolved_status. Usage: resolve_tokens.py"""
import pandas as pd
F = "/usr/local/jupyter/FVB_bias/revision_v2/freeze"
a = pd.read_csv(f"{F}/token_audit.tsv", sep="\t", dtype=str, keep_default_na=False)
r = pd.read_csv(f"{F}/recompute_report.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("claim_id")
MAP = {("LAR_main_text_v2.plain.txt", "1,047"): "K10", ("LAR_main_text_v2.plain.txt", "290"): "K11",
       ("response_to_reviewers_v2.plain.txt", "0.333"): "K16", ("response_to_reviewers_v2.plain.txt", "0.386"): "K17",
       ("cover_letter_v2.plain.txt", "0.333"): "K18", ("figures/Fig1_survey.pdf", "0.333"): "K12", ("figures/Fig1_survey.pdf", "0.386"): "K13",
       ("figures/Fig2_f1_ase.pdf", "257,514"): "G32", ("figures/Fig2_f1_ase.pdf", "32,067"): "G29", ("figures/Fig2_f1_ase.pdf", "0.483"): "G33",
       ("figures/SFig2_multimapping.pdf", "474"): "D6", ("LAR_main_text_v2.plain.txt", "±0.007"): "B4"}
a["resolved_by"] = [MAP.get((f, t), "") if c in ("UNMAPPED", "XLSX", "M2") else "" for f, t, c in zip(a.file, a.token, a["class"])]
a["resolved_status"] = [r.status.get(x, "") if x else "" for x in a.resolved_by]
# tokens whose value is reproduced by a PASS claim of the independent recomputation become class M2 (ref = claim id); tokens that cannot be
# recomputed (1,047 search hits; the 290-sample GEO series) stay UNMAPPED with their M-2 classification in resolved_status
m2 = (a.resolved_status == "PASS")
a.loc[m2, "ref"] = "M-2 claim " + a.loc[m2, "resolved_by"]; a.loc[m2, "class"] = "M2"
a.to_csv(f"{F}/token_audit.tsv", sep="\t", index=False)
u = a[a["class"] == "UNMAPPED"]
print(a["class"].value_counts().to_string()); print("UNMAPPED without an M-2 claim:", int((u.resolved_by == "").sum()))
print(u[["file", "token", "resolved_by", "resolved_status"]].to_string())
