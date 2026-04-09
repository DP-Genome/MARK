# Mitochondrial DNA Diagnostic Split-Track Pipelines

## Overview

This repository contains two mitochondrial DNA analysis pipelines built for control/hypervariable region using PowerSeq® CRM Nested System, Custom kit. It features a **diagnostic split-track design**:

- `ONT_MITO_v36.sh` for Oxford Nanopore Technologies data
- `ILM_MITO_v36.sh` for Illumina paired-end data

Both pipelines process targeted mtDNA control-region data against a **linearized mitochondrial reference** and then split each sample into two downstream analysis tracks:

1. **Baseline track**  
   No amplicon-boundary trimming is applied after alignment.
2. **Trimmed track**  
   Reads are reassigned to a single amplicon using a **maximum-overlap rule**, then trimmed in alignment space by rewriting the CIGAR with soft-clipping so that only the tiled target interval for that amplicon is retained.

Version 36 is a refinement of the recent split-track model and specifically introduces a fix for the **1 bp gap at the left border of amplicons** caused by 0-based coordinate conversion, while retaining the soft-clipping and maximum-overlap framework from recent versions.

---

## What changed in v36

### Shared v36 changes

- fixes the 1 bp gap on the left border of amplicons
- retains the v34/v35 soft-clipping strategy
- retains **maximum-overlap amplicon assignment** for read-to-amplicon selection
- uses **perfectly tiled amplicon bounds** for the trimmed track to eliminate internal gaps and overlaps between retained target intervals
- keeps the same dual-track comparison design so baseline and trimmed outputs can be compared directly

### Important terminology change

Older repository versions described the second track as a **masked** track. In v36, the scripts generate a **trimmed** track instead:

- output names now use `trimmed`
- BAM editing is performed by **CIGAR rewriting and soft-clipping**, not by simple base-quality masking terminology

---

## General pipeline design

Both pipelines follow the same high-level structure:

1. input discovery
2. initial read-retention safety filtering
3. QC and adapter/primer trimming
4. alignment to a linearized mtDNA reference
5. generation of two comparison tracks from the aligned BAM
6. independent variant calling on each track
7. annotation and export of multiple diagnostic VCFs
8. depth comparison between baseline and trimmed tracks

This repository is therefore designed not only to produce final call files, but also to make it easier to see **what changed, where reads were lost, and how the trimming logic affected depth and variants**.

---

## Input handling

### ONT pipeline

`ONT_MITO_v36.sh` accepts a folder containing:

- `.fastq`
- `.fastq.gz`

Each FASTQ is treated as one sample.

### Illumina pipeline

`ILM_MITO_v36.sh` accepts a folder containing paired-end FASTQs. The script searches for matching read pairs using names that contain either:

- `_R1` / `_R2`
- `.R1` / `.R2`

Each discovered pair is treated as one sample.

---

## Preprocessing and trimming

## 1. Length-safety filtering

Both pipelines begin with a length-safety step to remove clearly unsuitable reads before core analysis.

### ONT

Uses `fastplong`.

Default:

- `MAX_LEN=1500`

### Illumina

Uses `fastp` in paired mode before merging.

Default:

- `MAX_LEN=1500`

---

## 2. Platform-specific preprocessing

### ONT

After length filtering, the ONT pipeline optionally filters reads by the `qs:f:` tag in the FASTQ header.

Default:

- `QS_MIN=10`

If the `qs:f:` tag is absent, the pipeline keeps the reads and continues.

### Illumina

The Illumina pipeline merges paired reads with `fastp` and applies merged-read filtering during that step.

Default merge/filter settings:

- `READQ=20`
- `UNQUAL_PCT=40`
- `MIN_LEN=90`
- `N_BASE_LIMIT=5`

Because `--include_unmerged` is enabled, the merged FASTQ can contain both merged and retained unmerged reads, which is relevant when interpreting downstream read-retention counts.

---

## 3. QC

FastQC is run before alignment:

- ONT: on the QS-filtered FASTQ
- Illumina: on the merged FASTQ

A later FastQC run is also performed on the final **trimmed BAM** for diagnostic review.

---

## 4. Adapter and primer trimming

Both pipelines use a two-pass `cutadapt` strategy:

1. 3′ trimming with `-a file:<adapter_file>`
2. 5′ trimming with `-g file:<adapter_file>`

Optional fixed trimming from both ends is then applied with `EXTRA_TRIM`.

Defaults:

- `CUTADAPT_ERR=0.10`
- `CUTADAPT_OVL=5`
- `EXTRA_TRIM=0`
- `MIN_LEN=90`

After this, both pipelines perform an additional post-trim length filter using `fastp`.

Defaults:

- `MIN_LEN_POST=90`
- `MAX_LEN_POST=300`

---

## Reference and annotation inputs

Both pipelines expect:

- a linearized mtDNA reference FASTA
- a BED file describing annotated regions such as homopolymer regions and blacklist sites
- an adapter/primer list for `cutadapt`

Default filenames:

- `linearized_mtdna.fasta`
- `linearized_regions.bed`
- `Updated_Adapter_Primer_List_Cutadapt_cleaned.txt`

### Adapter file lookup behavior

If `ADAPTER_FILE` is not exported manually, the pipelines look for the adapter file in:

1. the parent directory of the input folder
2. the current working directory

### Automatic indexing

- the ONT pipeline auto-builds the minimap2 `.mmi` index if needed
- the Illumina pipeline auto-builds the `bwa-mem2` index if needed
- both pipelines create a FASTA index with `samtools faidx` if missing

---

## Alignment

### ONT

Alignment uses `minimap2` with the `map-ont` preset:

```bash
minimap2 -ax map-ont -t <threads> <ref.mmi> <reads.fastq>
```

### Illumina

Alignment uses `bwa-mem2`:

```bash
bwa-mem2 mem -t <threads> <ref.fasta> <reads.fastq>
```

The aligned reads are converted to BAM, sorted, and indexed before track splitting.

---

## Diagnostic split-track design

After the initial aligned BAM is created, each sample is split into two tracks.

### Track A: Baseline

The initial sorted BAM is copied directly:

- `*_baseline_sorted.bam`

This track serves as the unmodified comparison track.

### Track B: Trimmed

The initial BAM is streamed through an embedded Python routine that:

- skips unmapped and secondary alignments
- computes the reference span of each alignment from the CIGAR
- assigns each read to the amplicon with the **largest overlap**
- uses a fixed set of **tiled amplicon bounds** for retained intervals
- rewrites the alignment CIGAR with soft-clipping so only the assigned tiled interval is retained
- discards reads that fall completely outside the retained target interval
- strips leading and trailing deletions from the active retained segment
- writes the resulting alignments into:
  - `*_trimmed_sorted.bam`

This trimmed track is the experimental track used to examine the effect of amplicon-boundary enforcement on downstream calling.

---

## Amplicon logic used in v36

The embedded Python defines original amplicon spans and separate **tiled bounds** used for retention.

Original amplicon spans:

- Amp1: 7729 to 7842
- Amp2: 7832 to 7941
- Amp3: 7939 to 8124
- Amp4: 8103 to 8202
- Amp5: 8190 to 8315
- Amp6: 8271 to 8437
- Amp7: 8421 to 8542
- Amp8: 8531 to 8649
- Amp9: 8627 to 8721
- Amp10: 8714 to 8877

Tiled retention bounds used by the trimmed track:

- Amp1: 7729 to 7837
- Amp2: 7838 to 7940
- Amp3: 7941 to 8113
- Amp4: 8114 to 8196
- Amp5: 8197 to 8293
- Amp6: 8294 to 8429
- Amp7: 8430 to 8536
- Amp8: 8537 to 8638
- Amp9: 8639 to 8717
- Amp10: 8718 to 8877

These tiled bounds are what enforce the no-gap, no-overlap retained intervals in the trimmed track.

---

## Variant calling

Variant calling is performed independently on the baseline and trimmed BAMs.

Common logic:

```bash
bcftools mpileup -a FORMAT/AD,FORMAT/DP \
  -d <PILEUP_MAX_DEPTH> \
  -Q<BASEQ_MIN> \
  -q<MAPQ_MIN> \
  -Ou -f <ref> <bam> \
| bcftools call -mv --ploidy 1 -Ov
```

Defaults:

- `PILEUP_MAX_DEPTH=100000`
- `BASEQ_MIN=20`
- `MAPQ_MIN=20`
- ploidy `1`

Initial post-calling filter:

- `QUAL > SAFETY_QUAL`
- `INFO/DP >= MIN_DEPTH`

Defaults:

- `SAFETY_QUAL=20`
- `MIN_DEPTH=10`

---

## Annotation and VCF output generation

If the quality-filtered VCF is non-empty, it is annotated with the BED file so that variants can be labeled by region class.

Region labels currently documented by the header inserted during runtime:

- `HP_Region`
- `Blacklist_Site`

After annotation, each track produces several downstream VCFs.

### SNP-only VCF

Generated by retaining:

- `TYPE="snp"`
- `QUAL > STRICT_QUAL`

and then excluding the predefined artifact positions:

- `7898`
- `7899`
- `8595`

Outputs:

- `*_baseline_snps.vcf`
- `*_trimmed_snps.vcf`

### Clean VCF

Generated from SNPs passing the strict QUAL filter and then excluding:

- `RegionType="HP_Region"`
- `RegionType="Blacklist_Site"`

Outputs:

- `*_baseline_clean.vcf`
- `*_trimmed_clean.vcf`

### Homopolymer-only VCF

Retains variants annotated as:

- `RegionType="HP_Region"`

Outputs:

- `*_baseline_homopolymers.vcf`
- `*_trimmed_homopolymers.vcf`

### Additional intermediate VCFs

Each track also retains:

- `*_raw.vcf`
- `*_qual_filtered.vcf`
- `*_annotated_all.vcf`

These are diagnostic outputs and are useful when tracing why a variant appears or disappears between steps.

---

## Output files

For each sample, the most important output files are usually:

- `*_trimmed_snps.vcf` for the trimmed-track SNP callset
- `*_trimmed_clean.vcf` for the trimmed-track SNP callset after homopolymer and blacklist exclusion
- `*_baseline_snps.vcf` for direct comparison to the untrimmed baseline track

Because the repository is diagnostic by design, there is no single output that makes the others redundant. The files serve different purposes:

- `raw` and `qual_filtered` help inspect pre-annotation and early post-call filtering
- `annotated_all` shows region labeling
- `snps` gives the strict SNP-only view with artifact-position exclusion
- `clean` gives the stricter review set excluding homopolymer and blacklist regions
- `homopolymers` isolates calls in homopolymer regions for manual review

---

## Read-retention and depth diagnostics

Each run produces a run-level summary file:

- `run_summary.txt`

This records, per sample and per stage:

- reads in
- reads out
- dropped reads
- dropped percentage

Each sample also gets a log file:

- `<sample>.log`

This records command execution and includes a depth comparison summary:

- baseline average depth
- trimmed average depth

This makes v36 useful not just for generating callsets, but also for investigating where coverage is being lost.

---

## Output structure

Each run creates a timestamped output directory.

### ONT

```text
ONT_MITO_v36_<input_folder_name>_<timestamp>_output/
```

### Illumina

```text
ILM_MITO_v36_<input_folder_name>_<timestamp>_output/
```

Within that folder, each sample gets its own subdirectory.

Example:

```text
ONT_MITO_v36_FASTQS_20260408_120000_output/
└── barcode07/
    ├── barcode07_initial_sorted.bam
    ├── barcode07_baseline_sorted.bam
    ├── barcode07_trimmed_sorted.bam
    ├── barcode07_baseline_raw.vcf
    ├── barcode07_baseline_qual_filtered.vcf
    ├── barcode07_baseline_annotated_all.vcf
    ├── barcode07_baseline_snps.vcf
    ├── barcode07_baseline_clean.vcf
    ├── barcode07_baseline_homopolymers.vcf
    ├── barcode07_trimmed_raw.vcf
    ├── barcode07_trimmed_qual_filtered.vcf
    ├── barcode07_trimmed_annotated_all.vcf
    ├── barcode07_trimmed_snps.vcf
    ├── barcode07_trimmed_clean.vcf
    ├── barcode07_trimmed_homopolymers.vcf
    ├── barcode07.log
    └── ...
```

---

## Usage

### ONT

```bash
./ONT_MITO_v36.sh <input_fastq_folder>
```

### Illumina

```bash
./ILM_MITO_v36.sh <input_fastq_folder>
```

---

## Required software

### Common

- bash
- python3
- awk
- samtools
- bcftools
- cutadapt
- fastqc

### ONT

- minimap2
- fastplong

### Illumina

- bwa-mem2
- fastp

---

## Environment-variable configuration

The pipelines are designed to be configured through exported environment variables before execution.

### Common parameters

| Variable | Default | Description |
| --- | --- | --- |
| `threads` | `8` | Number of threads |
| `MIN_LEN` | `90` | Minimum length allowed during cutadapt trimming |
| `MIN_LEN_POST` | `90` | Minimum length after post-trim length filtering |
| `MAX_LEN_POST` | `300` | Maximum length after post-trim length filtering |
| `EXTRA_TRIM` | `0` | Fixed trimming from both read ends after cutadapt |
| `CUTADAPT_ERR` | `0.10` | Cutadapt error rate |
| `CUTADAPT_OVL` | `5` | Minimum overlap for cutadapt matching |
| `ref` | `linearized_mtdna.fasta` | Reference FASTA |
| `regions_bed` | `linearized_regions.bed` | BED file used for annotation |
| `ADAPTER_FILE` | auto-detected | Adapter/primer list for cutadapt |
| `SAFETY_QUAL` | `20` | Initial VCF QUAL threshold |
| `STRICT_QUAL` | `60` | Strict SNP QUAL threshold |
| `MIN_DEPTH` | `10` | Minimum depth required after calling |
| `PILEUP_MAX_DEPTH` | `100000` | Maximum mpileup depth |
| `BASEQ_MIN` | `20` | Minimum base quality for mpileup |
| `MAPQ_MIN` | `20` | Minimum mapping quality for mpileup |
| `EDGE_PROTECT` | `1` | Defined in the scripts, but not actively applied by the current v36 trimming logic |
| `DISCARD_WARN_PCT` | `5` | Defined threshold for warning use in future or external wrappers |

### ONT-specific parameters

| Variable | Default | Description |
| --- | --- | --- |
| `MAX_LEN` | `1500` | Maximum read length kept in the lensafe step |
| `QS_MIN` | `10` | Minimum `qs:f:` score if present |
| `mmi_index` | `${ref}.mmi` | Minimap2 index path |

### Illumina-specific parameters

| Variable | Default | Description |
| --- | --- | --- |
| `MAX_LEN` | `1500` | Maximum read length retained before merge |
| `READQ` | `20` | fastp quality threshold for merged-read filtering |
| `UNQUAL_PCT` | `40` | fastp unqualified percent limit |
| `N_BASE_LIMIT` | `5` | Maximum allowed number of `N` bases |

---

## Example usage

### ONT with defaults

```bash
./ONT_MITO_v36.sh /path/to/ont_fastqs
```

### ONT with custom QS and extra trimming

```bash
export QS_MIN=12
export EXTRA_TRIM=5
./ONT_MITO_v36.sh /path/to/ont_fastqs
```

### Illumina with custom merge settings

```bash
export READQ=25
export UNQUAL_PCT=30
export EXTRA_TRIM=5
./ILM_MITO_v36.sh /path/to/illumina_fastqs
```

### Custom reference and BED

```bash
export ref=/path/to/linearized_mtdna.fasta
export regions_bed=/path/to/linearized_regions.bed
export ADAPTER_FILE=/path/to/Updated_Adapter_Primer_List_Cutadapt_cleaned.txt
./ONT_MITO_v36.sh /path/to/ont_fastqs
```

---

## Interpretation notes

These pipelines are intentionally verbose in their outputs.

That is not bloat. It is the point of the design.

The large number of BAM and VCF outputs exists so you can answer practical questions such as:

- Did the variant disappear before or after annotation?
- Did it disappear only in the trimmed track?
- Was it removed because it was in a homopolymer region?
- Was it removed by the strict SNP filter?
- Did the trimmed track materially reduce average depth?
- Did a read-retention step remove more data than expected?

For that reason, the extra outputs should be viewed as:

- **comparison outputs** for baseline versus trimmed behavior
- **diagnostic outputs** for debugging read loss and variant loss
- **review outputs** for separating strict review callsets from region-flagged callsets

---

## Important notes

1. These pipelines use **linearized coordinates**. If downstream interpretation needs original mtDNA positions, a separate post-processing coordinate-reversion step is still required.
2. The variable `EDGE_PROTECT` is still present in both scripts, but the active v36 Python trimming logic uses direct tiled bounds and does **not** currently apply a dynamic edge-protection window.
3. The v36 scripts name the second track `trimmed`, not `masked`. Older README text that refers to `masked_*` outputs is no longer accurate for the current version.

---

## Summary

Version 36 documents a **dual-track, diagnostic mitochondrial DNA workflow** for ONT and Illumina data in which:

- preprocessing is platform-specific
- alignment is performed against a linearized mtDNA reference
- each sample is split into a baseline and a trimmed track
- the trimmed track uses maximum-overlap amplicon assignment and CIGAR soft-clipping against tiled amplicon intervals
- both tracks are called, filtered, annotated, and retained for side-by-side review

This makes the repository useful both for producing review-ready SNP callsets and for debugging how trimming and amplicon-boundary enforcement affect depth and variant retention.
