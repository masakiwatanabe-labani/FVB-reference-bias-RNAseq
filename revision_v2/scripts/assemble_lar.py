#!/usr/bin/env python3
"""Phase H-3: assemble the LAR/BMC-format main text and Additional file 1 from the section files, then build the docx files with
pandoc (--reference-doc = the v1 submission docx, so headings/body styles match).
  main text : title page (v1) -> structured Abstract (abstract_structured.md) -> Keywords -> Background (introduction.md, heading
              renamed) -> Methods -> Results -> Discussion -> Conclusions (split from discussion.md) -> List of abbreviations
              (abbreviations.md) -> Declarations (declarations.md) -> References (references.md, numbered by number_citations.py)
              -> Figure legends (main figures only) -> Table 1 (table1.md; footnote without script names)
  Additional file 1: Supplementary figure legends (S1–S3), Supplementary table captions, Supplementary Notes 1–2
Transformations applied to every section: local paths -> "the repository (file)"; numbered-list steps -> "Step n." paragraphs;
"* p < 0.05; ** p < 0.01; *** p < 0.001" asterisks escaped; HTML comments dropped; Markdown '---' rules dropped.
Outputs: revision_v2/submission/LAR_main_text_v2.md/.docx, LAR_additional_file_1.md/.docx (+ plain-text dumps for checks).
Usage: assemble_lar.py"""
import re, subprocess, os
PROJ="/usr/local/jupyter/FVB_bias"; M=f"{PROJ}/manuscript/revised_v2"; OUT=f"{PROJ}/revision_v2/submission"; REFDOC=f"{PROJ}/LAR text submit.docx"
os.makedirs(OUT, exist_ok=True)
def rd(f): return open(f"{M}/{f}", encoding="utf-8").read()
def clean(s):
    s=re.sub(r"<!--.*?-->", "", s, flags=re.S)
    s=re.sub(r"\n-{3,}\n", "\n", s)
    # local paths -> repository file names
    s=re.sub(r"`(?:revision_v2|analysis_v2|revision)/(?:[\w.\-]+/)*([\w.\-]+)`", r"the repository file `\1`", s)
    s=re.sub(r"\((?:revision_v2|analysis_v2)/(?:[\w.\-]+/)*([\w.\-]+)\)", r"(repository file \1)", s)
    s=s.replace("the repository file `manifest.tsv`","the sample manifest in the repository (`manifest.tsv`)")
    # numbered procedure lists -> Step paragraphs (Markdown "1. **X**: ..." at line start)
    s=re.sub(r"^(\d+)\.\s+\*\*([^*]+)\*\*:\s*", lambda m: f"Step {m.group(1)}. **{m.group(2)}**: ", s, flags=re.M)
    s=re.sub(r"^(\d+)\.\s+\*\*([^*]+)\*\*\s", lambda m: f"Step {m.group(1)}. **{m.group(2)}** ", s, flags=re.M)
    s=re.sub(r"\n(?!\n)(Step \d+\. )", r"\n\n\1", s)                      # each Step is its own paragraph (no soft line breaks)
    s=re.sub(r"`(?:revision_v2|analysis_v2)/(?:[\w.\-]+/)*`", "the repository", s)   # directory-only mentions
    s=re.sub(r"(?<![`\w/])(?:revision_v2|analysis_v2)/(?:[\w.\-]+/)*([\w.\-]+)?", lambda m: f"the repository file {m.group(1)}" if m.group(1) else "the repository", s)   # bare paths
    # significance asterisks (escape so pandoc does not treat them as emphasis)
    s=s.replace("* p < 0.05; ** p < 0.01; *** p < 0.001", "\\* p < 0.05; \\*\\* p < 0.01; \\*\\*\\* p < 0.001")
    return s.strip()+"\n"
def heading_shift(s, top):  # make the section's own H1 the given top-level heading text
    return re.sub(r"^# .*$", f"# {top}", s, count=1, flags=re.M)
title=open(f"{OUT}/title_page.md", encoding="utf-8").read()
abstract=rd("abstract_structured.md"); keywords=rd("keywords.md")
background=heading_shift(rd("introduction.md"), "Background"); methods=rd("methods.md"); results=rd("results.md")
disc=rd("discussion.md"); i=disc.find("## Conclusions"); discussion=disc[:i]; conclusions="# Conclusions\n"+disc[i+len("## Conclusions"):]
abbrev=rd("abbreviations.md"); decl=rd("declarations.md"); refs=rd("references.md")
legends=rd("figure_legends.md"); j=legends.find("**Supplementary Figure S1"); main_legends=heading_shift(legends[:j], "Figure legends"); supp_legends="# Supplementary figure legends\n\n"+legends[j:]
table1=rd("table1.md"); table1=re.sub(r"^Generated from the sample manifest.*\n", "", table1, flags=re.M)
main="\n\n".join(clean(x) for x in (title, abstract, keywords, background, methods, results, discussion, conclusions, abbrev, decl, refs, main_legends, table1))
paper_title=[l.strip() for l in title.splitlines() if l.strip() and not l.startswith("**") and not l.startswith("#")][0]   # first plain line of the title page = the title
add="\n\n".join(clean(x) for x in ("# Additional file 1\n\nSupplementary figure legends, supplementary table captions and Supplementary Notes 1–2 for: "+paper_title+"\n", supp_legends, rd("supplementary_table_captions.md"), rd("supplementary_note1.md"), rd("supplementary_note2.md")))
open(f"{OUT}/LAR_main_text_v2.md","w",encoding="utf-8").write(main); open(f"{OUT}/LAR_additional_file_1.md","w",encoding="utf-8").write(add)
for name in ("LAR_main_text_v2","LAR_additional_file_1"):
    subprocess.run(["pandoc",f"{OUT}/{name}.md","-f","markdown","-t","docx","--reference-doc",REFDOC,"-o",f"{OUT}/{name}.docx"],check=True)
    subprocess.run(["pandoc",f"{OUT}/{name}.docx","-t","plain","--wrap=none","-o",f"{OUT}/{name}.plain.txt"],check=True)
    print(name, "words:", len(open(f"{OUT}/{name}.md",encoding="utf-8").read().split()))
