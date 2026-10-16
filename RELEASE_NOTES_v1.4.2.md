# MARK v1.4.2

**Sample reports in forensic notation. The variant calls and VCFs are identical to v1.4.1.**

## What changed
The sample report (`Final_Pipeline_Results/sample_reports/<library>_final_sample_report.tsv`) writes variants the way
forensic mtDNA profiles are written (ISFG and EMPOP conventions):
- **SNVs:** `152C` (was `m.152C`). The HGVS form, `m.152T>C`, stays in the Notes column.
- **Insertions:** moved to the 3-prime end of the repeat they lie in and numbered there. The VCF record 310 T>TC
  becomes `315.1C`, and a three-C insertion in the 568–573 tract becomes `573.1C 573.2C 573.3C`.
- **Deletions:** also moved to the 3-prime end, one entry per deleted position. A CA deletion in the 514–523 repeat
  becomes `523del 524del`.
- **Traceability:** the Notes column keeps the VCF record of each insertion and deletion, e.g.
  `Insertion (VCF: 310 T>TC)`.

The VCFs are unchanged: they keep the standard left-aligned representation. Insertions and deletions remain
review-only (the "LENGTH VARIANTS AND HOMOPOLYMER SITES" section); the call set is still SNVs only.

## Known limitation
Around the HV2 C-stretch (303–315), freebayes can split one complex haplotype into two records, e.g. 310 T>C plus
an insertion. Each record is written on its own, which can give an unusual notation such as `310.1T 310.2C …`.
These sites are review-only; read them against the alignment, or align the haplotype phylogenetically (e.g. EMPOP).

## Validation
- **Unit test:** VCF records to notation for 310 T>TC → 315.1C, 152 T>C → 152C, +CCC at 567 → 573.1C–573.3C, GCA>G
  at 513 → 523del 524del, and C insertions in the HV1 C-tract → 16188.1C and 16193.1C.
- **Reference run:** reports rebuilt from the reference-run VCFs for 2800M (152C, 263G, 477C, 16519C; 315.1C under
  review), DNA 007 nanopore (573.1C–573.4C under review) and NIST C (523del 524del under review).
- **Full pipeline test** against the previous version, on DNA 007 · MinION 2025 SUP v5.2 (3 libraries), 2800M ·
  Illumina 2025 and NIST C · Illumina 2025: all 240 VCFs have identical records.
