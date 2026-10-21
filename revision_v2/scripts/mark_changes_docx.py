#!/usr/bin/env python3
"""Phase H-9: mark the text that changed relative to the submitted version in blue.
Compares revision_v2/submission/LAR_main_text_v2.docx (clean version, left untouched) with the v1 submission docx sentence by
sentence within each section (Heading 1 sections; Heading 2 subsections where present) — punctuation, whitespace, case and
citation numbers ([n], [n–m], [CITE …]) are ignored — and writes LAR_main_text_v2_marked.docx in which every sentence that
does not occur in v1 (new or changed) is blue (RGB 0,0,255). A unit (subsection, or section without subsections) in which
≥ REWRITE_FRAC of the sentences changed is marked blue as a whole with a one-line note at its top. Figure legends and Table 1
(table cells) are compared too; References are skipped (renumbered) with a note. A legend line is inserted at the top.
Usage: mark_changes_docx.py [v1.docx] [v2.docx] [out.docx] [--fuzzy R]
  defaults: 'LAR text submit.docx', LAR_main_text_v2.docx, *_marked.docx. Default rule = exact sentence identity after
  normalisation; --fuzzy R additionally treats a sentence as unchanged when its difflib similarity to a v1 sentence of the same
  section is ≥ R (e.g. 0.9: wording edits of a few words are not marked).
Re-run after every rebuild of the clean docx (build_submission.sh build calls it)."""
import re, sys, copy, difflib
FUZZY=float(sys.argv[sys.argv.index("--fuzzy")+1]) if "--fuzzy" in sys.argv else None
if "--fuzzy" in sys.argv: i=sys.argv.index("--fuzzy"); del sys.argv[i:i+2]
from docx import Document
from docx.shared import RGBColor
from docx.text.run import Run
PROJ="/usr/local/jupyter/FVB_bias"
V1=sys.argv[1] if len(sys.argv)>1 else f"{PROJ}/LAR text submit.docx"
V2=sys.argv[2] if len(sys.argv)>2 else f"{PROJ}/revision_v2/submission/LAR_main_text_v2.docx"
OUT=sys.argv[3] if len(sys.argv)>3 else V2.replace(".docx","_marked.docx")
BLUE=RGBColor(0,0,255); REWRITE_FRAC=0.6
H1=["Abstract","Keywords","Background","Introduction","Methods","Results","Discussion","Conclusions","List of abbreviations","Declarations","References","Figure legends","Table 1"]
def norm(s):
    s=re.sub(r"\[CITE[^\]]*\]","",s); s=re.sub(r"\[\d+(?:[–,\-]\d+)*\]","",s)
    return re.sub(r"[\W_]+","",s).lower()
ABBR=("Fig","Figs","al","vs","e.g","i.e","ca","no","approx","Supplementary Fig","Dr","St")
def sentences(text):
    """(start,end) spans of sentences in text"""
    spans=[]; start=0
    for m in re.finditer(r"[.!?;](?=\s+\S|$)", text):
        before=text[:m.start()]
        if any(before.endswith(a) for a in ABBR) or re.search(r"\d$", before) and text[m.start()]=="." and m.end()<len(text) and text[m.end():m.end()+2].strip()[:1].isdigit(): continue
        end=m.end(); spans.append((start,end)); start=end
    if start<len(text): spans.append((start,len(text)))
    return [(a,b) for a,b in spans if text[a:b].strip()]
# ---- v1: sentences per section (headings identified by exact text) --------------------------------------------------------------
v1=Document(V1); v1_sec={}; v1_all=set(); cur="(front)"
def h1_name(t):
    t=t.strip().rstrip(":"); return t if t in H1 else None
for p in v1.paragraphs:
    t=p.text.strip()
    if not t: continue
    n=h1_name(t)
    if n: cur=n; continue
    for a,b in sentences(t):
        k=norm(t[a:b])
        if len(k)>=6: v1_sec.setdefault(cur,set()).add(k); v1_all.add(k)
v1_cells=set()
for tb in v1.tables:
    for row in tb.rows:
        for c in row.cells:
            k=norm(c.text)
            if k: v1_cells.add(k)
# v1 Table 1 is plain text after the "Table 1." caption (no w:tbl): take every tab/newline-separated token as a cell
pars=v1.paragraphs
for i,p in enumerate(pars):
    if p.text.strip().startswith("Table 1"):
        for q in pars[i:i+40]:
            if q.text.strip().startswith("Supplementary"): break
            for tok in re.split(r"[\t\n]+", q.text):
                k=norm(tok)
                if k: v1_cells.add(k)
        break
print(f"v1: {len(v1_all)} sentences in sections {sorted(v1_sec)}; {len(v1_cells)} table cells")
# ---- v2: walk paragraphs, decide per sentence ----------------------------------------------------------------------------------
v2=Document(V2); units=[]; cur1="(front)"; cur2=None; unit=None
def new_unit(name, first_par):
    u={"name":name,"pars":[],"first":first_par,"n":0,"changed":0}; units.append(u); return u
for p in v2.paragraphs:
    st=p.style.name; t=p.text.strip()
    if st=="Heading 1" or (st=="Heading 2" and cur1 in ("(front)",)):
        cur1=t; cur2=None; unit=new_unit(cur1,p); continue
    if st=="Heading 2":
        cur2=t; unit=new_unit(f"{cur1} / {cur2}",p); continue
    if unit is None: unit=new_unit("(front)",p)
    if not t: continue
    if cur1=="References": continue
    pool=v1_sec.get(cur1,set())
    spans=[]; n=0; ch=0
    for a,b in sentences(p.text):
        k=norm(p.text[a:b])
        if len(k)<6: continue
        n+=1
        if k in pool or k in v1_all: continue
        if FUZZY and difflib.get_close_matches(k, list(pool) if pool else list(v1_all), n=1, cutoff=FUZZY): continue
        ch+=1; spans.append((a,b))
    unit["pars"].append((p,spans)); unit["n"]+=n; unit["changed"]+=ch
def color_runs(par, mask):
    """colour the characters of par where mask[i] is True, splitting runs at colour boundaries"""
    pos=0
    for r in list(par.runs):
        t=r.text; seg=mask[pos:pos+len(t)]; pos+=len(t)
        if not t: continue
        if all(seg): r.font.color.rgb=BLUE; continue
        if not any(seg): continue
        # split into maximal uniform pieces
        pieces=[]; i=0
        while i<len(t):
            j=i
            while j<len(t) and seg[j]==seg[i]: j+=1
            pieces.append((t[i:j],seg[i])); i=j
        r.text=pieces[0][0]
        if pieces[0][1]: r.font.color.rgb=BLUE
        prev=r._r
        for txt,col in pieces[1:]:
            nr=copy.deepcopy(r._r); prev.addnext(nr); run=Run(nr,par); run.text=txt
            run.font.color.rgb=BLUE if col else None
            prev=nr
def mark_all(par):
    for r in par.runs: r.font.color.rgb=BLUE
report=[]
for u in units:
    if u["name"]=="References":
        note=u["first"].insert_paragraph_before("[References are not marked: they were renumbered in order of first appearance / 参考文献は番号振り直しのため対象外]")
        # insert after the heading instead of before it
        u["first"]._p.addnext(note._p); mark_all(note); continue
    frac=u["changed"]/u["n"] if u["n"] else 0
    if u["n"]>0 and frac>=REWRITE_FRAC and u["pars"]:
        for p,_ in u["pars"]: mark_all(p)
        note=u["pars"][0][0].insert_paragraph_before("[This unit was rewritten in full / この節は全面的に書き直した]"); mark_all(note)
        report.append((u["name"],u["n"],u["changed"],"REWRITTEN"))
    else:
        for p,spans in u["pars"]:
            if not spans: continue
            mask=[False]*len(p.text)
            for a,b in spans:
                for i in range(a,b): mask[i]=True
            color_runs(p,mask)
        report.append((u["name"],u["n"],u["changed"],"sentence-level" if u["changed"] else "unchanged"))
# tables (Table 1): cells not present in v1 tables -> blue
tc=0; tn=0
for tb in v2.tables:
    for row in tb.rows:
        for c in row.cells:
            k=norm(c.text)
            if not k: continue
            tn+=1
            if k not in v1_cells:
                tc+=1
                for p in c.paragraphs: mark_all(p)
# legend at the very top
first=v2.paragraphs[0]; leg=first.insert_paragraph_before("Blue text = changed or added relative to the submitted version (青字＝提出版からの変更・追加). Comparison is sentence-wise within each section, ignoring punctuation, spacing and citation numbers; units rewritten in full carry a note at their top; References are not marked (renumbered)."); mark_all(leg)
v2.save(OUT)
tot_n=sum(r[1] for r in report); tot_c=sum(r[2] for r in report)
print(f"{OUT} (rule: {'fuzzy ≥ '+str(FUZZY) if FUZZY else 'exact sentence identity'})\nsentences compared {tot_n}, marked {tot_c} ({tot_c/tot_n:.0%}); table cells {tn}, marked {tc}")
for name,n,c,how in report: print(f"  {how:15s} {c:4d}/{n:<4d} {name}")
