# MARK v1.2.3 — changes since v1.2.1

This release also includes the changes made in v1.2.2. The variant-calling settings and filters are unchanged.

## Illumina: merged reads only (`MARK-I.sh`)

```
v1.2.1   merged reads (single-end)  +  unmerged pairs trimmed with cutadapt -a/-A, -g/-G and aligned as pairs
v1.2.3   merged reads (single-end) only
```

v1.2.0 began aligning the pairs fastp cannot merge. Part of that set is longer PCR product that runs from the forward primer of one amplicon to the reverse primer of a later one. At 300–550 bp it is too long to merge, and it appears in IGV as large-insert (red) pairs. In v1.2.3 only merged reads are aligned, so every aligned read is one collapsed fragment.

- fastp still writes unmerged pairs to `<sample>_unmerged_R1.fastq` / `_R2.fastq`, so they can be counted and inspected. The log reports how many there were.
- Reads whose mate failed quality filtering (`--unpaired1/--unpaired2`) are still kept as single-end reads, as in v1.2.1.
- The run summary no longer has a `2b_TrimPairs` row. `1_MergeFilter(Reads)` now counts merged reads only.
- Depth is lower where most fragments do not merge, typically across length-variable homopolymers. Check per-position depth there before relying on a call.

## Tool versions

The pipeline has been validated with:

| tool | version |
|---|---|
| fastp | 1.3.7 |
| fastplong | 0.7.1 |
| cutadapt | 5.2 |
| bwa-mem2 | 2.3 |
| minimap2 | 2.31 |
| samtools | 1.24 |
| bcftools | 1.24 |
| FastQC | 0.12.1 |

The conda recipe now sets these as minimum versions.

On macOS, cutadapt in multi-core mode can occasionally stop writing output part way through a large nanopore file, whatever the cutadapt version. The existing watchdog (`CUTADAPT_STALL_SECS`, default 600 s) catches this and reruns the step single-core. The output is identical; the step just takes longer.

## Linux: tool steps run in a plain shell (v1.2.2, `MARK.sh`, `MARK-I.sh`)

```
v1.2.1   run_log bash -lc "{ ... }"
v1.2.2   run_log bash -c  "{ ... }"
```

A login shell re-reads `/etc/profile`. On Debian and Ubuntu that file sets `PATH` outright, which drops the activated conda environment, so the alignment and calling steps could not find `bwa-mem2`, `samtools` or `bcftools`. On macOS it moved `/usr/local/bin` ahead of the environment, so a system copy of a tool could be used instead of the environment's. The pipelines no longer start login shells.
