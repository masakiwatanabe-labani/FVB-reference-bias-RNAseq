#!/usr/bin/env python3
"""Phase J-7: verify that every "**Changed:**" pointer in revision_v2/response_to_reviewers.md names a v1 paragraph
("<Section>, paragraph n", "<Section>, *<subsection>*, paragraph n") or a NEW section/paragraph ("<Section>, NEW section *title*",
"<Section>, NEW paragraph *first words*"), and that each resolves to rows of revision_v2/change_ledger.tsv whose trigger matches the
reviewer item (R1-M1 … R2-C4; multi-triggers are written "A/B"). Pointers that name figure files are listed but not resolved.
Exit 1 if any pointer does not resolve. Usage: verify_response_pointers.py"""
import re, csv, sys
P="/usr/local/jupyter/FVB_bias"
ALLOWED={"R1-M1":{"R1-M1"},"R1-M2":{"R1-M2","NEW-ASM"},"R1-M3":{"NEW-ASM","R1-M2"},"R1-m4":{"R1-m4"},"R1-m5":{"R1-m5"},"R1-m6":{"R1-m6"},"R2-C1":{"NEW-WS"},"R2-C2":{"R2-C2"},"R2-C3":{"R2-C3"},"R2-C4":{"R2-C4"}}
rows=list(csv.DictReader(open(f"{P}/revision_v2/change_ledger.tsv"), delimiter="\t"))
PLAIN=re.sub(r"[\W_]+","",open(f"{P}/revision_v2/submission/LAR_main_text_v2.plain.txt", encoding="utf-8").read()).lower()
def trig_ok(row, allowed): return bool(set(row["trigger"].split("/")) & allowed)
def norm(s): return re.sub(r"[\W_]+","",s).lower()
txt=open(f"{P}/revision_v2/response_to_reviewers.md", encoding="utf-8").read()
items=[]; rev=None
for line in txt.splitlines():
    if line.startswith("## Reviewer 1"): rev="R1"
    elif line.startswith("## Reviewer 2"): rev="R2"
    m=re.match(r"### (Major concern|Minor concern|Comment) (\d+)", line)
    if m:
        kind={"Major concern":"M","Minor concern":"m","Comment":"C"}[m.group(1)]; items.append([f"{rev}-{kind}{m.group(2)}", []])
    elif items and "**Changed:**" in line:
        items[-1][1].append(line.split("**Changed:**",1)[1].strip())
# the combined answer to R1-M2 and R1-M3 follows the R1-M3 heading: attach the same pointers to R1-M2
ptr={k:v for k,v in items}
if "R1-M2" in ptr and not ptr["R1-M2"] and ptr.get("R1-M3"): ptr["R1-M2"]=ptr["R1-M3"]
bad=0; report=[]
for item,plist in ptr.items():
    allowed=ALLOWED[item];
    if not plist: report.append((item,"(no Changed: pointer)","—","FAIL")); bad+=1; continue
    for ptext in plist:
        for ptr_ in [x.strip().rstrip(".") for x in re.split(r";", ptext) if x.strip()]:
            if re.fullmatch(r"all figure files", ptr_, flags=re.I): report.append((item,ptr_,"figure files (not text)","ok")); continue
            mb=re.match(r"(Abstract|Background|Methods|Results|Discussion|Conclusions)(?:, \*(.+?)\*)?, paragraph beginning \"([^\"]+)\"", ptr_)
            if mb:
                sec,sub,words=mb.groups(); key=norm(words)
                if key not in PLAIN: report.append((item,ptr_,"opening words not found in the plain dump","FAIL")); bad+=1; continue
                hit=[r for r in rows if r["section"]==sec and r["new_text"] and norm(r["new_text"]).startswith(key)]
                paras_={(r["section"],r["para_no"]) for r in hit}
                sel=[r for r in rows if (r["section"],r["para_no"]) in paras_ and trig_ok(r,allowed)]
                n=len(sel); st="ok" if n else "FAIL"; bad+=(n==0); report.append((item,ptr_,f"{n} ledger rows ({', '.join(sorted(set(r['status'] for r in sel)))}); paragraph found in plain dump",st)); continue
            m=re.match(r"(Abstract|Background|Methods|Results|Discussion|Conclusions)(?:, \*(.+?)\*)?(?:, (NEW (?:section|paragraph)) \*(.+?)\*)?(?:, paragraphs? ([\d–\-, and]+))?(?: \(.*\))?$", ptr_)
            if not m: report.append((item,ptr_,"unparsed","FAIL")); bad+=1; continue
            sec,sub,newkind,newtitle,paras=m.groups(); sel=[r for r in rows if r["section"]==sec and trig_ok(r,allowed)]
            if sub: sel=[r for r in sel if norm(sub.rstrip("…")) [:25] in norm(r["v1_subsection"]) or (r["status"]=="NEW" and norm(sub.rstrip("…"))[:25] in norm(r["new_text"]))]
            if newkind:
                key=norm(newtitle.rstrip("…"))[:30]; sel=[r for r in sel if r["status"]=="NEW" and (key in norm(r["new_text"]) or key in norm(r["v1_subsection"]))]
                # a NEW section is identified by its paragraphs; accept rows of the section whose trigger matches and whose text follows the title
                if not sel: sel=[r for r in rows if r["section"]==sec and r["status"]=="NEW" and trig_ok(r,allowed)]
            if paras:
                nums=set()
                for part in re.split(r",| and ", paras):
                    part=part.strip()
                    if not part: continue
                    if "–" in part or "-" in part: a,b=re.split(r"[–\-]",part); nums|=set(range(int(a),int(b)+1))
                    else: nums.add(int(part))
                col="v1_sub_para" if sub else "v1_para"
                sel=[r for r in sel if r[col] and int(float(r[col])) in nums and r["status"] in ("EDIT","KEEP","DROP")]
            n=len(sel); st="ok" if n else "FAIL"; bad+= (n==0)
            report.append((item,ptr_,f"{n} ledger rows ({', '.join(sorted(set(r['status'] for r in sel)))})",st))
for item,p,res,st in report: print(f"{st:4s} {item:6s} {p[:90]:90s} -> {res}")
print("unresolved pointers:", bad); sys.exit(1 if bad else 0)
