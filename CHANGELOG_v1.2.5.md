# MARK v1.2.5 — primer-pair check warns and never stops the run

This release includes everything in v1.2.4 (see `CHANGELOG_v1.2.4.md`). Nothing changes in trimming,
alignment, calling or filtering: the same reads reach the BAM as in v1.2.4. v1.2.4 remains available
as a fallback.

## What changed (TRIM_MODE=exact only)

In exact mode, reads without a primer pair of one amplicon are removed, exactly as before. What
changes is what happens when a library has only a small share of such reads.

| | v1.2.4 | v1.2.5 |
|---|---|---|
| reads without a primer pair | removed | removed (unchanged) |
| share below `PRIMER_MIN_PCT` (20%) | the whole run stops | setting removed: the run never stops |
| share below `PRIMER_WARN_PCT` | one line in the log | flagged prominently (below) |
| `PRIMER_WARN_PCT` default | 25% nanopore, 40% Illumina | **10% nanopore, 25% Illumina** |
| rejected reads | classified on a sample, then deleted | the same sample is kept, labelled, in a review FASTQ for **every** library |

Why: the share of reads with a primer pair varies widely between good libraries, because most of
what is sequenced is PCR by-product. Across the 89 validation libraries it ran from 6.2% to 55.9%
(median 16.5%), and a 20% stop would have halted 52 of 85 nanopore libraries. A stop also ended the
whole run, not just the library. The kept reads of a low-share library are clean, primer-verified
products, so there is no reason to discard them.

### How a flagged library is shown

- A banner in the terminal and in the library's log.
- `WARNINGS.txt` at the top of the run folder, one line per flagged library: its share, the flag
  level, the reads kept and the review file.
- The run summary row `2c_PrimerPairShare` becomes `2c_PrimerPairShare_LOW_see_WARNINGS`.
- A second banner at the end of the run lists every flagged library.

### Review files (every library)

Up to `PRIMER_DIAG_READS` (default 20,000) rejected reads per rejecting step, as they entered the
primer step, with the reason in the read header (`mark_reject=...`):

- Nanopore: `<library>_primer_rejects_review.fastq`, holding reads rejected by
  - the primer step: no primer found, only one primer of one amplicon, or primers from two amplicons
  - the adapter check: non-adapter sequence outside the primers, or adapter inside the product
- Illumina:
  - `<library>_primer_rejects_review_merged.fastq`: merged reads that are not one whole product
  - `<library>_primer_rejects_review_R1.fastq` and `_R2.fastq`: unmerged pairs whose mates do not
    start with the two primers of one amplicon

At the defaults each file is about 7–13 MB (measured on the test libraries). Illumina's pair check now writes rejected pairs to
a file instead of discarding them (`--untrimmed-output` in place of `--discard-untrimmed`). The
pairs that are kept are the same.

## Fixed

- **Reject reasons in the run summary.** In v1.2.4, rejected reads with no primer match at all were
  never counted, so "no primer found" was missing from `<library>_amplicon_trim.tsv`. Those reads
  were spread across "only one primer" and "primers from two amplicons". They are now counted
  (in one low-yield library, 1,680 of 33,363 rejected reads). Which reads are kept is unaffected.

## Verified

On a nanopore library, a low-yield flow-cell library and an Illumina library, v1.2.5 keeps exactly the
same reads as v1.2.4 (identical read names in the final BAM), and the snps and clean VCFs are identical.

- The low-yield library (6.2% of reads carrying a primer pair) was flagged and still processed.
- The Illumina library (43%) was flagged only when the flag level was raised to 50% for the test.

## Also

- Version 1.2.5 in `MARK.sh`, `MARK-I.sh`, `MARKLaunch.py` and `conda-recipe/meta.yaml`.
