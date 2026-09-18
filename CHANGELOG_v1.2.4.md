# MARK v1.2.4 — exact amplicon trimming

This release also includes the changes made in v1.2.2 and v1.2.3 (see `CHANGELOG_v1.2.3.md`). The
variant-calling settings and filters are unchanged.

## A second trimming mode: `TRIM_MODE`

```
TRIM_MODE=adapters   the validated path, trimming with the adapter lists   (default)
TRIM_MODE=primers    cut every read at its own amplicon primers
```

Trimming with adapter lists can never be exact by construction. An adapter entry that stops one base
short of the junction leaves a base on the read; an entry that contains primer bases cuts into the
amplicon. In `primers` mode a read is kept only if it carries the forward **and** reverse primer of the
**same** amplicon, with nothing but library adapter outside them, and is cut exactly at the outer
primer ends. The mapped span of a full-length product is then the product size in
`CRM_Nested_products.bed` by construction.

- The cutadapt primer files are written at run time from `CRM_Nested_primers_empirical.bed` and the
  reference, so they cannot drift from the coordinates the rest of the pipeline uses.
- The primers stay on the read (`--action=retain`) on purpose: the product coordinates include them,
  and the existing primer-footprint step soft-clips them in the BAM, so no primer base reaches the
  caller.
- Each read name gains the amplicon its primer pair identified (`..._Amp7`), so amplicon assignment no
  longer depends on alignment coordinates alone.
- The run summary reports, per amplicon, reads kept and reads rejected with the reason (no primer pair,
  only one primer, primers from two amplicons, non-adapter sequence outside the primers, adapter inside
  the product).

### Illumina (`MARK-I.sh`), primers mode

1. fastp merge as before, but **no mate rescue**: a fragment is kept only if both mates survive.
2. Pairs that did not merge are quality-trimmed at the 3' end (`--cut_right` Q20) and merged again.
   (`--cut_right` in the first pass is worse on every amplicon, so it is used only here.)
3. Merged reads must be one whole product: `cutadapt -g file:<anchored> --revcomp --action=retain`,
   where each entry is `^FORWARD...REVERSE$`. Anchoring is what rejects read-through products; without
   it a product spanning two amplicons would be trimmed down into a perfect-looking shorter one.
4. Pairs that still did not merge are kept only if read 1 and read 2 start with the two primers of one
   amplicon (`--pair-adapters`), and are aligned as a pair.
5. The `-a`/`-g` adapter-list passes are not used in this mode.

### Nanopore (`MARK.sh`), primers mode

The `-b` adapter-list pass is replaced by the linked-primer pass (`-e 0.15 -O 18`), followed by a filter
that requires the sequence outside each primer to be the library adapter junction, and rejects a product
that contains adapter **inside** it (a chimera whose junction falls between two primer matches).
`minimap2` gains `--end-bonus 15`: with exact trimming the read ends are real template, and without it
the aligner soft-clips error-containing primer ends (exact spans 88.9% → 96.6%, calls unchanged).

### Guard

If fewer than `PRIMER_MIN_PCT` (default 20) per cent of reads carry a matching primer pair the run
aborts, so a kit or primer change fails visibly instead of silently discarding everything; below
`PRIMER_WARN_PCT` (default 25) it warns and continues. On the validation libraries the share is 24–43%,
because most of what is sequenced is by-product; two dimer-heavy nanopore libraries sit at 16–17%.

## Validation

All 89 libraries, both platforms, against the validated v1.2.1 runs:

- Illumina: reported calls identical in 4/4 libraries; allele fractions move by at most 0.001.
- Nanopore: identical in 79 of 81 real libraries. The two exceptions are flongle libraries with mean
  depth 72x and 20x that had already failed.
- Amplicon sizes match the BED for all 10 amplicons on both platforms; exact spans 98.7–98.9%
  (Illumina) and 96.3% (nanopore).
- No library adapter survives trimming, and no read-through product reaches the BAM.
- Cost: nanopore keeps a median 55% of reads and 72% of mean depth; the loss is by-product, and
  concentrated in the amplicon overlap windows.
