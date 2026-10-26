#!/usr/bin/env python3
"""Phase W-3 (2026-10-26): Table 1 as a stand-alone docx fragment for the authors to transplant into their manuscript (the manuscript
is maintained outside the repository, so assemble_lar.py is not run). Input: manuscript/revised_v2/generated/table1.md (make_table1.py)
with the "Generated from ..." line dropped, as assemble_lar.py does; pandoc with the same reference docx; landscape page; the table
layout of assemble_lar.py (fix_tables, taken from its source without executing the assembly).
Output: manuscript/revised_v2/generated/table1_fragment.docx (+ .plain.txt). Usage: make_table1_docx.py"""
import ast, re, subprocess
PROJ = "/usr/local/jupyter/FVB_bias"; G = f"{PROJ}/manuscript/revised_v2/generated"; REFDOC = f"{PROJ}/LAR text submit.docx"
OUT = f"{G}/table1_fragment.docx"
src = open(f"{PROJ}/revision_v2/scripts/assemble_lar.py", encoding="utf-8").read()
fn = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == "fix_tables")
ns = {}; exec(compile(ast.Module(body=[fn], type_ignores=[]), "assemble_lar.py:fix_tables", "exec"), ns)
md = re.sub(r"^Generated from the sample manifest.*\n", "", open(f"{G}/table1.md", encoding="utf-8").read(), flags=re.M)
open(f"{G}/table1_fragment.md", "w", encoding="utf-8").write(md)
subprocess.run(["pandoc", f"{G}/table1_fragment.md", "-f", "markdown", "-t", "docx", "--reference-doc", REFDOC, "-o", OUT], check=True)
from docx import Document
from docx.enum.section import WD_ORIENT
d = Document(OUT); sec = d.sections[-1]
if sec.orientation != WD_ORIENT.LANDSCAPE:
    sec.orientation = WD_ORIENT.LANDSCAPE; sec.page_width, sec.page_height = sec.page_height, sec.page_width
d.save(OUT)
print("tables:", ns["fix_tables"](OUT))
subprocess.run(["pandoc", OUT, "-t", "plain", "--wrap=none", "-o", OUT.replace(".docx", ".plain.txt")], check=True)
