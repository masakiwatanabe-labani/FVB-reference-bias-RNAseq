"""Phase Z: reference-genome classification of a series from the full !Sample_data_processing text of its RNA-Seq samples
(replaces the keyword classify() of stage3_survey_v2.py, which missed strain-aware methods written in other words and mis-read some
assembly names). Case-insensitive regular expressions, applied per sentence:
  strain_aware (priority): STRAIN_AWARE patterns, in a sentence about alignment, mapping, quantification, index or reference
                           (CONTEXT), so that a sentence naming only the strain of the animals does not count;
  stated_B6:               B6 patterns (GRCh / hg are human and never count);
  not_reported:            the rest; reference_named_without_assembly = Y when an Ensembl/RefSeq/UCSC/iGenomes release or a
                           GRCh name is given without a mouse assembly name.
classify_reference(lines) -> dict(category, matched_keyword, sentence, assemblies_named, strain_aware_method,
                                  reference_named_without_assembly)"""
import re
F = re.I
STRAIN_AWARE = [
    (r"FVB[_/ \-]?N[J]?[_\-]{0,2}(v\d|R\b|R\d|genome|transcriptom|cDNA|reference|assembl)", "FVB/NJ genome or transcriptome"),
    (r"Mus_musculus_Fvb", "FVB/NJ genome or transcriptome"),
    (r"GCA_001624535", "FVB/NJ genome or transcriptome"),   # Phase Z2: GenBank accession of FVB_NJ_v1 (stated as such in GSE135432)
    (r"pseudo-?genome", "strain-specific pseudogenome"),
    (r"personali[sz]ed", "personalised genome"),
    (r"individuali[sz]ed genome", "personalised genome"),
    (r"diploid genome", "personalised genome"),
    (r"strain[- ]specific (genome|reference|transcriptome|pseudogenome)", "strain-specific pseudogenome"),
    (r"\bWASP\b", "allele-aware filtering"),
    (r"g2gtools", "personalised genome"),
    (r"\bN-?mask", "allele-aware filtering"),
    (r"SNPsplit", "allele-aware filtering"),
    (r"Seqnature", "strain-specific pseudogenome"),
    (r"\bEMASE\b", "allele-aware filtering"),
    (r"allele[- ]specific (align|mapping)", "allele-aware filtering"),
    (r"strain[- ]aware", "other"),
    (r"SNP[- ]aware", "allele-aware filtering"),
    (r"FVB.{0,40}(genome|transcriptome|cDNA|reference)|(genome|transcriptome|cDNA|reference).{0,40}FVB", "FVB/NJ genome or transcriptome"),
]
CONTEXT = r"align|map|quantif|index|reference|genome|transcriptom|cDNA|assembl|build|salmon|kallisto|\bstar\b|hisat|bowtie|tophat|bwa|rsem|cell ?ranger|gsnap|pseudo"
B6 = [r"GRCm\s?\.?\s?3[89]", r"GRCme8", r"\bmm\s?(9|10|39)\b", r"GENCODE\s?(release\s?)?v?M\d+|\bRelease M\d+", r"C57BL/?6", r"\bB6 (genome|reference)",
      r"NCBIM?3[78]", r"NCBI\s?(build\s?)?m?3[78]\b", r"MGSCv?3[67]",
      # Phase Z2
      r"\bmm8\b", r"(mouse|mus musculus|\bmm)[^.]{0,40}\bbuild\s?3[78]\b", r"MOUSE_B3[78]", r"GRCm38\.\d{2}", r"ENSEMBL_38_\d{2}", r"GENCODE[^.]{0,20}\bM\d+\b|\bM\d+ (genome )?annotation"]
MOUSE_SENT = r"mouse|mus musculus|murine|\bmm\b"
ENSEMBL_SENT = r"ensembl?e?"
AMBIGUOUS = r"\bmm\s?38\b|GRChm38"
NO_ASSEMBLY = r"ensembl|refseq|ucsc|igenomes?|GRCh\s?3[78]|\bhg(19|38)\b"
ASSEMBLY = [("GRCm38", r"GRCm\s?\.?\s?38|GRCme8|\bmm\s?10\b|NCBIM?38|NCBI\s?(build\s?)?m?38\b|(mouse|mus musculus|\bmm)[^.]{0,40}\bbuild\s?38\b|MOUSE_B38|ENSEMBL_38_\d{2}"),
            ("GRCm39", r"GRCm\s?\.?\s?39|\bmm\s?39\b"),
            ("mm9", r"\bmm\s?9\b|NCBIM?37|NCBI\s?(build\s?)?m?37\b|MGSCv?37|(mouse|mus musculus|\bmm)[^.]{0,40}\bbuild\s?37\b|MOUSE_B37"), ("mm8", r"\bmm8\b")]
HUMAN = r"GRCh\s?3[78]|\bhg\s?(19|38)\b"
def sentences(lines):
    out = []
    for l in lines:
        out += [s.strip() for s in re.split(r"(?<=[.!?;])\s+(?=[A-Z(])", l) if s.strip()]
    return out
def tag(s):
    """Phase Z2: mark tokens that name a mouse build only in context: B37/B38 in a sentence naming mouse -> MOUSE_B37/38;
    "38.NN" in a sentence naming Ensembl and mouse -> ENSEMBL_38_NN (Ensembl's GRCm38.NN release naming)"""
    if re.search(MOUSE_SENT, s, F):
        s = re.sub(r"\bB(3[78])\b", r"MOUSE_B\1", s)
        if re.search(ENSEMBL_SENT, s, F): s = re.sub(r"(?<![\w.])38\.(\d{2})\b", r"ENSEMBL_38_\1", s)
    return s
def method_of(sent, label):
    t = sent.lower()
    if "g2gtools" in t or re.search(r"personali[sz]ed|individuali[sz]ed genome|diploid genome", t): return "personalised genome (g2gtools)" if "g2gtools" in t else "personalised genome"
    if re.search(r"pseudo-?genome|seqnature|strain[- ]specific (genome|reference|pseudogenome)", t): return "strain-specific pseudogenome"
    if re.search(r"\bwasp\b|snpsplit|\bemase\b|n-?mask|allele[- ]specific (align|mapping)|snp[- ]aware", t): return "allele-aware filtering (WASP/SNPsplit)"
    if "gca_001624535" in t and not re.search(r"transcriptom|cdna|salmon|kallisto", t): return "FVB/NJ genome assembly"
    if re.search(r"fvb", t):
        return "FVB/NJ transcriptome or cDNA" if re.search(r"transcriptom|cdna|salmon|kallisto", t) else "FVB/NJ genome assembly"
    return "other"
def gencode_assembly(t):
    out = set()
    for m in re.finditer(r"(?:GENCODE[^.]{0,20}\b|GENCODE\s?(?:release\s?)?v?|\bRelease )M(\d+)\b|\bM(\d+)(?= (?:genome )?annotation)", t, F):
        n = int(m.group(1) or m.group(2))
        out.add("GRCm38" if n <= 25 else "GRCm39")
    return out
def _untag(res):
    res["matched_keyword"] = res["matched_keyword"].replace("MOUSE_B", "B").replace("ENSEMBL_38_", "38."); return res
def classify_reference(lines):
    return _untag(_classify(lines))
def _classify(lines):
    S0 = sentences(lines); S = [tag(x) for x in S0]; ORIG = dict(zip(S, S0)); res = dict(category="not_reported", matched_keyword="", sentence="", assemblies_named="", strain_aware_method="", reference_named_without_assembly="N", human_assembly_named="N", assembly_name_ambiguous="N")
    for s in S:
        if not re.search(CONTEXT, s, F): continue
        for pat, _ in STRAIN_AWARE:
            m = re.search(pat, s, F)
            if m:
                res.update(category="strain_aware", matched_keyword=m.group(0), sentence=ORIG[s], strain_aware_method=method_of(s, _)); break
        if res["category"] == "strain_aware": break
    asm = set()
    for s in S:
        for a, pat in ASSEMBLY:
            if re.search(pat, s, F): asm.add(a)
        asm |= gencode_assembly(s)
    res["assemblies_named"] = ";".join(a for a in ("GRCm38", "GRCm39", "mm9", "mm8") if a in asm)
    if res["category"] == "strain_aware": return res
    for s in S:
        for pat in B6:
            m = re.search(pat, s, F)
            if m: res.update(category="stated_B6", matched_keyword=m.group(0), sentence=ORIG[s]); return res
    for s in S:   # Phase Z2: a spelling that could be mm10/GRCm38 or an error is flagged, not decided
        m = re.search(AMBIGUOUS, s, F)
        if m: res.update(assembly_name_ambiguous="Y", matched_keyword=m.group(0), sentence=ORIG[s]); return res
    for s in S:
        m = re.search(NO_ASSEMBLY, s, F)
        if m and re.search(CONTEXT, s, F):
            res.update(reference_named_without_assembly="Y", matched_keyword=m.group(0), sentence=ORIG[s])
            if re.search(HUMAN, s, F): res["human_assembly_named"] = "Y"
            break
    return res
PATTERN_LISTS = {"strain_aware": [p for p, _ in STRAIN_AWARE], "strain_aware_context": CONTEXT, "stated_B6": B6, "assembly_normalisation": ASSEMBLY,
                 "reference_named_without_assembly": NO_ASSEMBLY, "human_assembly_named": HUMAN, "assembly_name_ambiguous": AMBIGUOUS, "context_tags (B37/B38 with mouse; 38.NN with Ensembl and mouse)": [MOUSE_SENT, ENSEMBL_SENT]}
