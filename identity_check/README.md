# Sample-identity check: the fixed 10% site subset

Every animal except the F1 hybrids (the animals of the three-arm analysis and the four rejected GSE55788 animals) was checked
against its recorded strain with GATK ASEReadCounter at a fixed random 10% subset of the identity-check sites; the F1 animals
were evaluated at the full set.

- `identity_sites_sub10pct.seed20261016.tsv` — the subset actually used (CHROM, POS, REF, ALT; GRCm39; ALT = FVB/NJ allele).
- `identity_intervals_sub10pct.seed20261016.bed` — the traversal intervals passed to ASEReadCounter with `-L` (subset sites
  merged within 1 kb); they only limit where the walker goes, the per-site counts are unchanged.
- `make_identity_subset.py` — regenerates both from the full site set; its output reproduces the subset VCF used in the analysis
  record for record and the interval file byte for byte. The full set (203,510 sites) is built by
  `revision_v2/scripts/A2_build_qc_sites.sh` followed by `bcftools +setGT … -t a -n c:0/1` (command in the script header).
- `subset_design_and_validation_2026-10-16.md` — the record written when the subset was introduced: definition
  (`random.Random(20261016).sample(range(203510), 20351)`, file order kept) and the validation against the full set (13 animals,
  identical verdicts).

Verdict rules: `revision_v2/scripts/strain_identity_qc.py` / `analysis_v2/pipeline/identity_verdict.py`; per-animal results:
`analysis_v2/results/qc/sample_identity_verification.tsv` (column `site_set`).
