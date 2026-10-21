#!/usr/bin/env python3
"""Markdown -> .docx with python-docx (no pandoc on this machine). Handles: #/##/### headings, paragraphs, **bold**, *italic*,
`code`, > block quotes, -/1. lists, | tables |, --- rules (page/section break), ``` code blocks, HTML comments (dropped),
$$ display equations (rendered as a centred text line: \\text{}, \\log_2, _{sub}, \\Delta, \\_ handled).
Usage: md_to_docx.py <in.md> <out.docx> [title]"""
import re, sys
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT

def eq_runs(par, s):
    """display equation -> runs with real subscripts (\\text{}, \\log_2, _{...}, bare _x, \\Delta, \\_)"""
    s=s.replace("\\log_2","log₂").replace("\\Delta","Δ").replace("\\_","\x02").replace("\\text{","\x01").replace(" - "," − ")
    while "\x01" in s: s=re.sub(r"\x01([^{}\x01]*)\}", r"\1", s)
    s=s.replace("\\","").replace("\x02","_"); s=re.sub(r"\s+"," ",s).strip()
    for tok in re.split(r"(_\{[^{}]*\}|_[A-Za-z]+)", s):
        if not tok: continue
        if tok.startswith("_"): r=par.add_run(tok[2:-1] if tok.startswith("_{") else tok[1:]); r.font.subscript=True; r.italic=True
        else: r=par.add_run(tok); r.italic=True

BOLD=re.compile(r"(\*\*.+?\*\*)"); INLINE=re.compile(r"(\*[^*]+\*|`[^`]+`)")
def add_runs(par, text, base_bold=False, base_italic=False):
    for seg in BOLD.split(text):
        if not seg: continue
        if seg.startswith("**") and seg.endswith("**") and len(seg)>4: add_runs(par, seg[2:-2], base_bold=True, base_italic=base_italic); continue
        for tok in INLINE.split(seg):
            if not tok: continue
            if tok.startswith("*") and tok.endswith("*") and len(tok)>2: r=par.add_run(tok[1:-1]); r.italic=True; r.bold=base_bold
            elif tok.startswith("`") and tok.endswith("`"): r=par.add_run(tok[1:-1]); r.font.name="Consolas"; r.font.size=Pt(9.5); r.bold=base_bold
            else: r=par.add_run(tok); r.bold=base_bold; r.italic=base_italic

def convert(src, dst, title=None):
    doc=Document(); st=doc.styles["Normal"]; st.font.name="Calibri"; st.font.size=Pt(11)
    for s in doc.sections: s.left_margin=s.right_margin=Cm(2.5); s.top_margin=s.bottom_margin=Cm(2.5)
    lines=open(src,encoding="utf-8").read().replace("\r","").split("\n")
    i=0; n=len(lines); para=[]; in_comment=False; in_code=False; code=[]
    def flush():
        nonlocal para
        if para:
            p=doc.add_paragraph(); add_runs(p," ".join(x.strip() for x in para)); p.paragraph_format.space_after=Pt(8); para=[]
    while i<n:
        l=lines[i]; s=l.strip()
        if in_comment:
            if "-->" in s: in_comment=False
            i+=1; continue
        if s.startswith("<!--"):
            flush()
            if "-->" not in s: in_comment=True
            i+=1; continue
        if s.startswith("```"):
            flush()
            if in_code: p=doc.add_paragraph(); r=p.add_run("\n".join(code)); r.font.name="Consolas"; r.font.size=Pt(9); code=[]; in_code=False
            else: in_code=True
            i+=1; continue
        if in_code: code.append(l); i+=1; continue
        if not s: flush(); i+=1; continue
        if s.startswith("$$"):
            flush(); eq=s.strip("$ ")
            if not eq:  # multi-line equation
                i+=1; buf=[]
                while i<n and not lines[i].strip().startswith("$$"): buf.append(lines[i]); i+=1
                eq=" ".join(buf)
            p=doc.add_paragraph(); eq_runs(p, eq); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; i+=1; continue
        m=re.match(r"^(#{1,4})\s+(.*)$", s)
        if m:
            flush(); lvl=len(m.group(1)); h=doc.add_heading("",level=min(lvl,4)); add_runs(h, m.group(2)); i+=1; continue
        if re.match(r"^-{3,}$", s): flush(); doc.add_paragraph().add_run("").add_break(); i+=1; continue
        if s.startswith("|"):
            flush(); rows=[]
            while i<n and lines[i].strip().startswith("|"):
                r=lines[i].strip().strip("|"); cells=[c.strip() for c in r.split("|")]
                if not all(re.match(r"^:?-{2,}:?$", c) for c in cells): rows.append(cells)
                i+=1
            if rows:
                ncol=max(len(r) for r in rows); t=doc.add_table(rows=len(rows), cols=ncol); t.style="Table Grid"; t.alignment=WD_TABLE_ALIGNMENT.CENTER
                for ri,r in enumerate(rows):
                    for ci in range(ncol):
                        cell=t.cell(ri,ci); cell.text=""; p=cell.paragraphs[0]; add_runs(p, r[ci] if ci<len(r) else "", base_bold=(ri==0))
                        for run in p.runs: run.font.size=Pt(8.5)
                doc.add_paragraph()
            continue
        if s.startswith(">"):
            flush(); buf=[]
            while i<n and lines[i].strip().startswith(">"): buf.append(lines[i].strip()[1:].strip()); i+=1
            p=doc.add_paragraph(); p.paragraph_format.left_indent=Cm(1); add_runs(p," ".join(buf), base_italic=True); continue
        m=re.match(r"^(\s*)([-*]|\d+\.)\s+(.*)$", l)
        if m:
            flush(); style="List Number" if m.group(2)[0].isdigit() else "List Bullet"
            p=doc.add_paragraph(style=style); add_runs(p, m.group(3)); i+=1
            while i<n and lines[i].startswith("   ") and lines[i].strip() and not re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]): add_runs(p," "+lines[i].strip()); i+=1
            continue
        para.append(s); i+=1
    flush(); doc.save(dst); print(f"{dst}: {len(doc.paragraphs)} paragraphs, {len(doc.tables)} tables")

if __name__=="__main__": convert(sys.argv[1], sys.argv[2])
