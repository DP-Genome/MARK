# MARK v1.6

**Heteroplasmy report.** Every allele at 5% or more of the reads is listed, with its fraction, in a VCF and in the
sample report. The calls are identical to v1.5.

## What changed
- **Why a second freebayes run:** the calls come from a one-copy model, which reports the majority base at each
  position. A 10% heteroplasmy at a non-variant site, or the minor allele of a 70/30 site, is never shown.
- **What v1.6 runs:** freebayes again on the final (primer-clipped) BAM, in pooled-continuous mode. This mode gives
  allele fractions, not genotypes.
- **What is reported:** every allele that meets all three of these:
  - fraction of at least `HET_MIN_AF`;
  - fraction of at most 1 − `HET_MIN_AF`;
  - depth of `MIN_DEPTH` or more.
- **Where it is reported:**
  - **`<library>_trimmed_heteroplasmy.vcf`:** the fraction is in `INFO/HET_AF`, and homopolymer regions and excluded
    sites are annotated (`RegionType`). After post-processing it is in `corrected_vcfs/` (rCRS).
  - **The sample report:** a new section, `### HETEROPLASMY: ALLELES AT 5% OR MORE (for review, not calls) ###`. Each
    entry is point (PHP) or length (LHP) heteroplasmy, in forensic notation, with homopolymer regions and excluded
    sites marked.
- **Complex alleles:** freebayes sometimes writes an insertion next to a base change as one allele. Its single-base
  parts are added up, so the fraction is the share of all reads at the position that carry the base.
- **New setting `HET_MIN_AF`:** default `0.05`; `0` turns the report off. It is in the launcher's parameter table
  (with help text), `--help`, the README settings table and the run summary.
- **freebayes requirement:** the report needs freebayes. With `VARIANT_CALLER=bcftools` and no freebayes installed,
  the report is skipped with a note in the log.

## Unchanged
Trimming, alignment, calls (`_snps`, `_clean`), filters, and the other outputs are the same as v1.5.

## Validation
Investigation: `analysis_and_comparisons/v1.5/investigations/heteroplasmy_indels_20261027`.
Test runs: `analysis_and_comparisons/v1.6/tests/heteroplasmy_20261027`.

**In-silico mixtures** of 2800M and DNA 007 reads (2–50%, both directions, Illumina and nanopore, 18 mixtures):

| Heteroplasmies | Reported |
|---|---|
| True fraction ≥ 5% (110) | 109 (the one missed was at 5.0%); fractions within 0.8 percentage points |
| True fraction 2% | none |

Mutserve 2.0.3 gives the same results, but only on alignments with reads on both strands. MARK's alignments have
every read on one strand, and on them Mutserve reports no heteroplasmy at all. LoFreq and Mutect2 are less accurate.

**Real libraries:**
- **2800M:** no heteroplasmy on either platform. The one exception is 16373C at 5–7%, which appears in every sample of
  the 2025 MinION run (a run artefact).
- **DNA 007:**
  - 385G at 5–8% in every library on both platforms. This is a probable low-level heteroplasmy; other samples are
    at ≤ 0.3%.
  - 460C at about 45% (Illumina) and about 70% (nanopore), from the C8 tract.
  - The 573 C-tract length heteroplasmy (+4C dominant).
  - Amp10 read-2 noise at 437–498 on Illumina.

**End-to-end runs:** command line, then `--postprocess`.
- **Nanopore:** 2800M and DNA 007, SUP v5.2 (2025), 6 libraries. Exit 0.
- **Illumina:** DNA007-R1 and 2800M. Exit 0.
- **Outputs:** heteroplasmy VCFs are in `corrected_vcfs/`, and the sample reports have the new section. The values
  match the standalone analysis.
