#!/usr/bin/env python3
"""Phase B-3: per-design counts, FASTQ volume, run-time and disk estimates from manifest.tsv.
Timing basis (old logs, 8 threads): naive STAR 1.30 min and FVB STAR 0.74 min per million read-pairs
(median over 18 fructose-study samples, see revision_v2/logs/B3 output); fastp + WASP filter + featureCounts
+ QC are taken as 0.6 min per million read-pairs (old stage5 logs). Disk per sample ≈ 2.4× FASTQ bytes
(trimmed FASTQ + 3 BAMs + indices; old data/ shows ~2.2–2.6×)."""
import csv, collections
PROJ="/usr/local/jupyter/FVB_bias"
M=list(csv.DictReader(open(f"{PROJ}/revision_v2/manifest/manifest.tsv"), delimiter="\t"))
def est(rows):
    reads=sum(int(r["read_count_total"]) for r in rows)/1e6; by=sum(int(r["fastq_bytes_total"]) for r in rows)
    mins=reads*(1.30+0.74+0.6)
    return len(rows), reads, by/1e9, mins/60, by*2.4/1e9
for name,flt in [("Design 1 (main): water only", lambda r:r["role"]=="strain_comparison" and r["design_main"]=="Y"),
                 ("Design 2 (sensitivity): water + fructose", lambda r:r["role"]=="strain_comparison" and r["design_sensitivity"]=="Y"),
                 ("Negative control (PyVT, pending)", lambda r:r["role"]=="negative_control_candidate")]:
    rows=[r for r in M if flt(r)]
    n,reads,gb,hrs,disk=est(rows)
    print(f"\n## {name}\n animals={n}  read-pairs={reads:.0f}M  FASTQ={gb:.1f} GB  3-arm wall time≈{hrs:.1f} h (8 threads, sequential)  disk≈{disk:.0f} GB")
    c=collections.Counter((r["tissue"],r["strain"],r["diet"]) for r in rows)
    for k,v in sorted(c.items()): print("  ",k,v)
