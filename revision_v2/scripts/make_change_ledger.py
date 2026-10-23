#!/usr/bin/env python3
"""Phase J — change ledger. Compares every sentence of the reconstructed sections (manuscript/revised_v2/*.md, written by
build_phase_j.py) with the submitted text (revision_v2/submission/LAR_v1_text.md; Supplementary Notes from manuscript/
supplementary_note{1,2}.md) and writes revision_v2/change_ledger.tsv:
  section, para_no, sent_no, status (KEEP/EDIT/NEW/DROP), trigger, v1_text, new_text, v1_section, v1_subsection, v1_para, number_check
KEEP = identical after normalisation (punctuation, whitespace, case and citation numbers ignored); EDIT = closest v1 sentence of the
same section has difflib similarity ≥ 0.6 and the paragraph trigger is not a NEW-* trigger; NEW otherwise; DROP = v1 sentence that is
neither kept nor the source of an edit. Triggers come from revision_v2/ledger_paragraph_triggers.tsv (one per paragraph).
number_check: numbers with ≥ 3 significant digits in EDIT/NEW sentences are looked up among the formatted values of numbers.tsv;
unmatched numbers are listed for manual review. Usage: make_change_ledger.py"""
import re, csv, difflib, collections
P="/usr/local/jupyter/FVB_bias"; M=f"{P}/manuscript/revised_v2"
H1=["Abstract","Keywords","Introduction","Methods","Results","Discussion","Conclusions","List of abbreviations","Declarations","References","Figure legends","Table caption","Supplementary Table Captions"]
def norm(s):
    s=re.sub(r"\[CITE[^\]]*\]","",s); s=re.sub(r"\\\[\d+(?:[–,\-]\d+)*\\\]","",s); s=re.sub(r"\[\d+(?:[–,\-]\d+)*\]","",s)
    s=s.replace("\\","").replace("−","-").replace("–","-").replace("--","-")
    return re.sub(r"[\W_]+","",s).lower()
ABBR=("Fig","Figs","al","vs","e.g","i.e","ca","no","approx","Dr","St")
def sentences(text):
    spans=[]; start=0
    for m in re.finditer(r"[.!?;:](?=\s+\S|$)", text):
        before=text[:m.start()]
        if any(before.endswith(a) for a in ABBR): continue
        if text[m.start()]==":" and not text[m.end():].strip().startswith(("bias =","alt_count","Δ =","residual","$$","Naive","FVB","samtools","- **")): continue
        end=m.end(); spans.append((start,end)); start=end
    if start<len(text): spans.append((start,len(text)))
    return [text[a:b].strip() for a,b in spans if text[a:b].strip()]
# ---- v1 ------------------------------------------------------------------------------------------------------------------------
V1=[]  # dicts: section, subsection, para (within section), sub_para, text
sec="(front)"; sub=""; np_=0; nsp=0
for line in open(f"{P}/revision_v2/submission/LAR_v1_text.md", encoding="utf-8"):
    t=line.strip()
    if not t: continue
    m=re.fullmatch(r"\*\*(.+?)\*\*:?", t)
    if m and m.group(1) in H1: sec=m.group(1); sub=""; np_=0; nsp=0; continue
    if m or re.fullmatch(r"\*\*\*.+\*\*\*", t): sub=re.sub(r"\*","",t).rstrip(":"); nsp=0; continue
    if sec in ("References","(front)") : continue
    np_+=1; nsp+=1
    for s in sentences(t): V1.append({"section":"Background" if sec=="Introduction" else sec,"subsection":sub,"para":np_,"sub_para":nsp,"text":s,"norm":norm(s)})
for nn,f in ((1,"supplementary_note1.md"),(2,"supplementary_note2.md")):
    np_=0
    for para in re.split(r"\n\s*\n", open(f"{P}/manuscript/{f}", encoding="utf-8").read()):
        t=para.strip()
        if not t or t.startswith("#") or t.startswith("|") or t.startswith("```") or t.startswith("$$"): continue
        np_+=1
        for s in sentences(t): V1.append({"section":f"Supplementary Note {nn}","subsection":"","para":np_,"sub_para":np_,"text":s,"norm":norm(s)})
v1_by_norm=collections.defaultdict(list)
for i,r in enumerate(V1): v1_by_norm[r["norm"]].append(i)
# ---- numbers.tsv formatted values -----------------------------------------------------------------------------------------------
NUM=set()
for r in csv.DictReader(open(f"{M}/generated/numbers.tsv"), delimiter="\t"):
    try: x=float(r["value"])
    except: continue
    for f in (".0f",".1f",".2f",".3f",".4f"):
        NUM.add(format(x,f).lstrip("-+")); NUM.add(format(abs(x),f))
    NUM.add(f"{int(round(abs(x))):,}"); NUM.add(str(int(round(abs(x)))))
    NUM.add(f"{x:.1e}".split("e")[0]); NUM.add(f"{x/1e6:.1f}"); NUM.add(f"{x*100:.1f}"); NUM.add(f"{x*100:.2f}"); NUM.add(f"{x*100:.0f}")
    if 0<x<100: NUM.add(f"{x:.0f}")
IGNORE={"20261016","0.50","0.5","0.25","1.0","0.3","0.30","0.1","0.9","0.10","0.90","0.55","2021","2026","110","14","20","255","100","4","8","3","0.05","0.70","80","1,047","2112"}
def num_check(s):
    s=re.sub(r"GSE\d+|GSM\d+|SRR\d+|ENSMUSG\d+|v\d+(\.\d+)*|Fig\. ?S?\d+[A-D]?|Table S?\d+[a-f]?|Note \d|REL-\d+|GRCm3\d|mm\d+|chromosome \d+|10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+","",s)
    bad=[]
    for n in re.findall(r"(?<![A-Za-z\d_])\d[\d,]*\.?\d*(?![A-Za-z\d_])", s):
        n=n.rstrip(".").rstrip(",");
        if n in IGNORE or len(n.replace(",","").replace(".","").lstrip("0"))<3: continue
        if n not in NUM: bad.append(n)
    return ";".join(bad)
# ---- new text --------------------------------------------------------------------------------------------------------------------
TR={(r["section"],int(r["para_no"])):r["trigger"] for r in csv.DictReader(open(f"{P}/revision_v2/ledger_paragraph_triggers.tsv"), delimiter="\t")}
FILES=[("abstract_structured.md","Abstract"),("keywords.md","Keywords"),("introduction.md","Background"),("methods.md","Methods"),("results.md","Results"),("discussion.md","Discussion"),("figure_legends.md","Figure legends"),("supplementary_table_captions.md","Supplementary Table Captions"),("supplementary_note1.md","Supplementary Note 1"),("supplementary_note2.md","Supplementary Note 2"),("declarations.md","Declarations")]
rows=[]; used=set()
for fname,secname in FILES:
    paras=[x.strip() for x in re.split(r"\n\s*\n", open(f"{M}/{fname}", encoding="utf-8").read()) if x.strip()]
    cur=secname; generated=False
    for i,ptxt in enumerate(paras,1):
        trig=TR.get((fname,i),"?")
        if trig=="H":
            if fname=="discussion.md" and ptxt.startswith("## Conclusions"): cur="Conclusions"
            if ptxt.startswith("<!-- BEGIN"): generated=True
            if ptxt.startswith("<!-- END"): generated=False
            continue
        if ptxt.startswith("<!--"): continue
        for j,s in enumerate(sentences(ptxt),1):
            k=norm(s)
            if len(k)<4: continue
            same=[x for x in v1_by_norm.get(k,[]) if V1[x]["section"]==cur] or v1_by_norm.get(k,[])
            if same:
                x=same[0]; used.add(x); v=V1[x]
                rows.append([cur,i,j,"KEEP",trig,v["text"],s,v["section"],v["subsection"],v["para"],v["sub_para"],""]); continue
            status="NEW"; vtext=""; vs=vsub=vp=vsp=""
            if not trig.startswith("NEW-"):
                pool=[x for x in range(len(V1)) if V1[x]["section"]==cur] or list(range(len(V1)))
                best=max(pool, key=lambda x: difflib.SequenceMatcher(None,k,V1[x]["norm"]).ratio()) if pool else None
                if best is not None and difflib.SequenceMatcher(None,k,V1[best]["norm"]).ratio()>=0.6:
                    status="EDIT"; used.add(best); v=V1[best]; vtext=v["text"]; vs,vsub,vp,vsp=v["section"],v["subsection"],v["para"],v["sub_para"]
            rows.append([cur,i,j,status,trig,vtext,s,vs,vsub,vp,vsp,"" if generated else num_check(s)])
for x,v in enumerate(V1):
    if x in used or v["section"] in ("Keywords","List of abbreviations","Declarations"): continue
    rows.append([v["section"],"","","DROP","CORR",v["text"],"",v["section"],v["subsection"],v["para"],v["sub_para"],""])
# DROP triggers by v1 subsection / paragraph (default CORR)
DROPTRIG={("Results","***Fau* and *Rsph3a* as worked examples of reference-genome bias**"):"NEW-WE",("Discussion",6):"NEW-WE",("Discussion",7):"NEW-WE",("Discussion",8):"R2-C3",("Discussion",9):"R2-C3",("Discussion",10):"R2-C3",("Discussion",11):"R2-C3",("Discussion",17):"CORR",("Discussion",18):"NEW-NC",("Discussion",19):"NEW-NC"}
for r in rows:
    if r[3]=="DROP":
        r[4]=DROPTRIG.get((r[7],r[9]), DROPTRIG.get((r[7],r[8]), "CORR"))
        if r[7]=="Results" and ("Rsph3a" in r[5] or "Fau" in r[5] or "Atp5g2" in r[5]): r[4]="NEW-WE"
        if ("The Ensembl 110 GTF was provided with" in r[5] and "sjdbGTFfile" in r[5]) or r[5]=="WASP SAM-tag output was enabled during this alignment.": r[4]="CLARITY"
        if "PyVT" in r[5] or "pvt_discordant" in r[5] or "copy number" in r[5] or "copy-number" in r[5]: r[4]="NEW-NC"
        if r[7]=="Results" and r[5].startswith(("This near-zero genome-wide effect","Instead, it reflects the low density","In F1 hybrid analyses, where reads","In pure-strain comparisons, the proportion")): r[4]="R2-C3/PHASE-P"   # Phase P-4
with open(f"{P}/revision_v2/change_ledger.tsv","w",encoding="utf-8",newline="") as o:
    w=csv.writer(o, delimiter="\t"); w.writerow(["section","para_no","sent_no","status","trigger","v1_text","new_text","v1_section","v1_subsection","v1_para","v1_sub_para","number_check"]); w.writerows(rows)
cnt=collections.defaultdict(lambda: collections.Counter())
for r in rows: cnt[r[0]][r[3]]+=1
print(f"{'section':30s} KEEP  EDIT  NEW   DROP")
for s in ["Abstract","Keywords","Background","Methods","Results","Discussion","Conclusions","Figure legends","Table caption","Supplementary Table Captions","Supplementary Note 1","Supplementary Note 2","Declarations"]:
    c=cnt.get(s)
    if c: print(f"{s:30s} {c['KEEP']:4d}  {c['EDIT']:4d}  {c['NEW']:4d}  {c['DROP']:4d}")
tot=collections.Counter(r[3] for r in rows); print("total", dict(tot))
body=[r for r in rows if r[0] in ("Abstract","Background","Methods","Results","Discussion","Conclusions") and r[3]!="DROP"]
print(f"body sentences {len(body)}: KEEP {sum(r[3]=='KEEP' for r in body)} ({sum(r[3]=='KEEP' for r in body)/len(body):.0%})")
nc=[(r[0],r[11]) for r in rows if r[11]]; print("number_check flags:", len(nc));
for s,n in nc[:40]: print("  ", s, n)
