#!/usr/bin/env python3
"""
check_figure_order.py — extract the order of first citation of every figure and
panel in a Markdown manuscript, and list violations of numerical order.

Usage:
  python3 check_figure_order.py MANUSCRIPT.md [--all-sections] [--legend-check]

By default only the main text (everything before a line starting with
'# Figure Legends') is scanned; --all-sections scans the whole file.
--legend-check additionally parses the Figure Legends section for panel labels
(**(A)** ...) and reports panels that are never cited in the main text.

Recognised citation forms (case-insensitive): Fig. 3A | Figure 3A | Fig 3A |
Figs. 3A and 4B | Fig. 3A,B | Fig. 3A, B | Fig. 3A–C | Fig. 3B, C |
Supplementary Figure S1 | Supplementary Fig. S1A | Fig. S1A | Figure S1A,B
"""
import re, sys, argparse

# figure word is case-insensitive; panel letters must be UPPERCASE (so that
# "Fig. 4A, highlighted" does not yield a panel "h")
CITE_RE = re.compile(
    r"(?i:(?:Supplementary\s+)?Fig(?:ure|s?\.?))\s*([Ss]?\d+)\s*"
    r"([A-Z](?:\s*(?:,|–|-|\band\b)\s*[A-Z](?![a-z]))*)?(?![a-z])")
PANEL_SPLIT = re.compile(r"[A-Z]")


def expand_panels(s):
    """'A,B' -> [A,B]; 'A–C' -> [A,B,C]; 'A and B' -> [A,B]; 'B, C' -> [B,C]"""
    if not s:
        return []
    out = []
    # ranges first
    for m in re.finditer(r"([A-Z])\s*[–-]\s*([A-Z])", s):
        a, b = m.group(1), m.group(2)
        out.extend(chr(c) for c in range(ord(a), ord(b) + 1))
    s2 = re.sub(r"([A-Z])\s*[–-]\s*([A-Z])", " ", s)
    out.extend(PANEL_SPLIT.findall(re.sub(r"\band\b", " ", s2)))
    seen, res = set(), []
    for p in out:
        if p not in seen:
            seen.add(p); res.append(p)
    return res


def fig_sort_key(f):
    return (1 if f.startswith("S") else 0, int(f.lstrip("S")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manuscript")
    ap.add_argument("--all-sections", action="store_true")
    ap.add_argument("--legend-check", action="store_true")
    ap.add_argument("--stop-at", default="# Figure Legends",
                    help="line prefix at which main-text scanning stops")
    args = ap.parse_args()

    lines = open(args.manuscript, encoding="utf-8").read().splitlines()
    main_lines = []
    for i, l in enumerate(lines, 1):
        if not args.all_sections and l.strip().startswith(args.stop_at):
            break
        main_lines.append((i, l))

    first_fig, first_panel, order = {}, {}, []
    for i, l in main_lines:
        for m in CITE_RE.finditer(l):
            fig = m.group(1).upper()
            panels = expand_panels(m.group(2))
            if fig not in first_fig:
                first_fig[fig] = (i, m.group(0).strip())
                order.append(fig)
            for p in panels:
                key = (fig, p)
                if key not in first_panel:
                    first_panel[key] = (i, m.group(0).strip())

    print("== First citation of each figure (main text order) ==")
    for f in order:
        i, txt = first_fig[f]
        print(f"  Fig {f:>3}  line {i:4d}  '{txt}'")

    print("\n== First citation of each panel ==")
    for (f, p), (i, txt) in sorted(first_panel.items(), key=lambda kv: (fig_sort_key(kv[0][0]), kv[0][1])):
        print(f"  Fig {f}{p}  line {i:4d}  '{txt}'")

    violations = []
    # 1. figure-level order (main figures and supplementary figures checked separately)
    for prefix in ("", "S"):
        figs = [f for f in order if f.startswith("S") == (prefix == "S")]
        nums = [int(f.lstrip("S")) for f in figs]
        for k in range(1, len(nums)):
            if nums[k] < max(nums[:k]):
                violations.append(f"Figure {figs[k]} first cited (line {first_fig[figs[k]][0]}) after a higher-numbered figure ({prefix}{max(nums[:k])})")
    # 2. panel-level order within each figure
    for f in order:
        panels = sorted([(p, first_panel[(ff, p)][0]) for (ff, p) in first_panel if ff == f], key=lambda x: x[1])
        letters = [p for p, _ in panels]
        for k in range(1, len(letters)):
            if letters[k] < max(letters[:k]):
                violations.append(f"Fig {f}{letters[k]} first cited (line {panels[k][1]}) after Fig {f}{max(letters[:k])}")
        # missing panels between A and max cited letter
        if letters:
            for c in range(ord("A"), ord(max(letters)) + 1):
                if chr(c) not in letters:
                    violations.append(f"Fig {f}{chr(c)} is never cited in the main text (but Fig {f}{max(letters)} is)")

    if args.legend_check:
        legend_panels = {}
        cur = None
        for l in lines:
            m = re.match(r"\*\*(?:Supplementary\s+)?Figure\s+(S?\d+)\.", l)
            if m:
                cur = m.group(1).upper(); legend_panels.setdefault(cur, set()); continue
            if cur:
                for pm in re.finditer(r"\*\*\(([A-Z])\)\*\*", l):
                    legend_panels[cur].add(pm.group(1))
        for f, ps in sorted(legend_panels.items(), key=lambda kv: fig_sort_key(kv[0])):
            if f not in first_fig:
                violations.append(f"Figure {f} (in legends) is never cited in the main text")
            for p in sorted(ps):
                if (f, p) not in first_panel:
                    violations.append(f"Fig {f}{p} (in legend) is never cited in the main text")

    print("\n== Violations ==")
    if violations:
        for v in violations: print("  !!", v)
        print(f"\n{len(violations)} violation(s)")
        sys.exit(1)
    else:
        print("  none — all figures and panels first appear in numerical order")
        sys.exit(0)


if __name__ == "__main__":
    main()
