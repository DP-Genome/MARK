# Mitochondrial DNA Diagnostic Split-Track Pipelines v32

## Overview

This repository contains two specialized bioinformatics pipelines for mitochondrial DNA (mtDNA) analysis using a **diagnostic split-track design**:

- **`ONT_MITO_v32.sh`** for Oxford Nanopore Technologies (ONT) data
- **`ILM_MITO_v32.sh`** for Illumina paired-end data

These pipelines are designed for **targeted mitochondrial control region analysis** using a **linearized mtDNA reference**, with a particular focus on handling **overlapping amplicon edge effects** and evaluating the impact of read masking on downstream variant calling.

Unlike a standard single-track workflow, each pipeline generates **two parallel analysis tracks per sample**:

1. **Baseline track**  
   Standard alignment and variant calling without amplicon masking.

2. **Masked track**  
   Experimental masking logic is applied to aligned reads using midpoint-based amplicon assignment, with a minimal boundary rescue window controlled by `EDGE_PROTECT`.

The **main intended product** of these pipelines is:

- **`*_masked_snps.vcf`**

All other outputs are produced for **evaluation, troubleshooting, validation, and double-checking**, including comparison against the baseline track and inspection of filtering behavior.

Both pipelines also generate **run-level read-retention summaries** so that read loss can be tracked step by step through preprocessing, trimming, and filtering.

---

## Primary Output

The main output file for each sample is:

- **`*_masked_snps.vcf`**

This file represents the **masked-track SNP-only VCF** after:

- alignment
- masking
- variant calling
- QUAL and depth filtering
- annotation
- exclusion of specific artifact positions

This is the file intended to serve as the primary variant callset for downstream review and interpretation.

---

## General Pipeline Design

Both v32 pipelines retain the original midpoint-assignment masking model from version 31, but introduce a small boundary rescue window to reduce loss of true variants near assigned amplicon edges.

### Key v32 changes

- Retains original midpoint-based read-to-amplicon assignment
- Adds a minimal **`EDGE_PROTECT`** window at assigned amplicon boundaries
- Preserves the masking model elsewhere to avoid broad callset disruption
- Removes reference FASTA loading from Python memory to improve speed

---

## Pipeline Workflow

The general steps are summarized below.

---

## Step 1: Input Handling

### ONT pipeline
`ONT_MITO_v32.sh` accepts a folder containing `.fastq` or `.fastq.gz` files. Each file is treated as one sample and processed independently.

### Illumina pipeline
`ILM_MITO_v32.sh` accepts a folder containing paired-end FASTQ files. It searches for matching `R1` and `R2` pairs and processes each sample independently.

---

## Step 2: Read-Length Safety Filtering

Both pipelines apply an initial read-length control step before core processing.

### ONT
Reads longer than the configured `MAX_LEN` threshold are removed using `fastplong`.

- Default: `MAX_LEN=1500`

### Illumina
R1 and R2 reads are filtered with `fastp` length-limit mode before merging.

- Default: `MAX_LEN=3000`

This step is intended to remove anomalous or unsuitable reads before further analysis.

---

## Step 3: Read Quality and Preprocessing

### ONT
After length filtering, the ONT pipeline optionally filters reads based on the `qs:f:` tag in the FASTQ header.

- Default: `QS_MIN=10`

If the `qs:f:` tag is absent, the pipeline keeps the reads and proceeds without QS filtering.

### Illumina
The Illumina pipeline performs merged-read generation using `fastp` with:

- quality threshold filtering
- unqualified base percentage control
- minimum merged-read length filtering
- `N` base limit filtering

Default Illumina merge/filter parameters:

- `READQ=20`
- `UNQUAL_PCT=40`
- `MIN_LEN=90`
- `N_BASE_LIMIT=5`

---

## Step 4: Quality Control

FastQC is run on the processed FASTQ input prior to alignment.

- ONT: FastQC is run on the QS-filtered FASTQ
- Illumina: FastQC is run on the merged FASTQ

An additional FastQC run is also performed on the final masked BAM output for diagnostic review.

---

## Step 5: Adapter and Primer Trimming

Both pipelines use a two-pass `cutadapt` strategy:

1. **3′ trimming** using `-a file:<adapter_file>`
2. **5′ trimming** using `-g file:<adapter_file>`

Optional fixed-end trimming is then applied using `EXTRA_TRIM`.

- Default: `EXTRA_TRIM=0`
- Minimum retained read length after trimming: `MIN_LEN=90`

The adapter/primer/barcode file is expected to be:

`Updated_Adapter_Primer_List_Cutadapt_cleaned.txt`

The pipeline looks for this file in:
- the parent directory of the input folder
- the current working directory

Or it can be set manually with `ADAPTER_FILE`.

---

## Step 6: Alignment

### ONT
Alignment is performed with `minimap2` using the `map-ont` preset against the linearized mitochondrial reference.

```bash
minimap2 -ax map-ont reference.mmi sample.fastq | samtools sort -o sample_initial_sorted.bam
```

### Illumina
Alignment is performed with `bwa-mem2` against the linearized mitochondrial reference.

```bash
bwa-mem2 mem reference.fasta sample.fastq | samtools sort -o sample_initial_sorted.bam
```

Reference indexing is generated automatically if missing.

---

## Step 7: Diagnostic Split-Track Generation

After the initial sorted BAM is created, the pipeline splits into two tracks.

### Track A: Baseline
The initial BAM is copied directly and used as the unmasked comparison track.

Outputs from this track are used for:
- internal comparison
- diagnostic review
- confirmation that masking is not removing expected signal

### Track B: Masked
The initial BAM is streamed through an embedded Python masking step that:

- assigns each read to the nearest amplicon by midpoint
- preserves bases within the assigned amplicon
- applies a minimal rescue window at the boundaries using `EDGE_PROTECT`
- masks bases outside the retained region by lowering base qualities
- removes reads falling completely outside the retained assigned region

This masked BAM is then used for downstream calling.

The key purpose of this step is to reduce overlap-related artifacts while preserving genuine edge variants as much as possible.

---

## Step 8: Variant Calling

Variant calling is performed independently for both baseline and masked tracks using `bcftools mpileup` and `bcftools call`.

Common defaults:

- ploidy: `1`
- `FORMAT/AD` and `FORMAT/DP` retained
- `PILEUP_MAX_DEPTH=100000`
- `BASEQ_MIN=20`
- `MAPQ_MIN=20`

Initial post-calling filter:

- `QUAL > SAFETY_QUAL`
- `INFO/DP >= MIN_DEPTH`

Default values:

- `SAFETY_QUAL=20`
- `MIN_DEPTH=10`

---

## Step 9: Annotation and Output Splitting

Filtered VCFs are annotated using the supplied BED file to mark regions such as:

- `HP_Region`
- `Blacklist_Site`

Each track then produces multiple VCF outputs.

### Main output
- **`*_masked_snps.vcf`**

### Additional masked outputs
- `*_masked_annotated_all.vcf`
- `*_masked_clean.vcf`
- `*_masked_homopolymers.vcf`
- `*_masked_raw.vcf`
- `*_masked_qual_filtered.vcf`

### Baseline comparison outputs
- `*_baseline_snps.vcf`
- `*_baseline_annotated_all.vcf`
- `*_baseline_clean.vcf`
- `*_baseline_homopolymers.vcf`
- `*_baseline_raw.vcf`
- `*_baseline_qual_filtered.vcf`

### SNP filtering details
The SNP-only files are restricted to SNPs passing the stricter threshold:

- `QUAL > STRICT_QUAL`

Default:

- `STRICT_QUAL=60`

The SNP-only outputs also exclude specific predefined artifact positions:

- `7898`
- `7899`
- `8595`

### Clean VCF logic
The `*_clean.vcf` outputs exclude calls annotated as:

- `HP_Region`
- `Blacklist_Site`

### Homopolymer VCF logic
The `*_homopolymers.vcf` outputs retain only variants annotated as:

- `HP_Region`

---

## Diagnostic Outputs

The v32 pipelines are designed to support troubleshooting and validation, not only final calling.

### 1. Run summary
Each run produces:

- `run_summary.txt`

This file tracks read retention and read loss across preprocessing stages for every sample.

### 2. Sample log
Each sample folder contains:

- `*.log`

This log records the executed commands and pipeline progression.

### 3. Depth comparison
For each sample, the pipeline calculates and prints:

- baseline average depth
- masked average depth

This is intended to help assess the impact of the masking step on usable coverage.

---

## Output Structure

Each run produces a timestamped output directory:

### ONT
```text
ONT_MITO_v32_<input_folder_name>_<timestamp>_output/
```

### Illumina
```text
ILM_MITO_v32_<input_folder_name>_<timestamp>_output/
```

Within each run folder, each sample gets its own subdirectory containing intermediate and final outputs.

Example:

```text
ONT_MITO_v32_FASTQS_20260326_120000_output/
└── barcode07/
    ├── barcode07_masked_snps.vcf
    ├── barcode07_masked_clean.vcf
    ├── barcode07_masked_homopolymers.vcf
    ├── barcode07_masked_annotated_all.vcf
    ├── barcode07_baseline_snps.vcf
    ├── barcode07_baseline_clean.vcf
    ├── barcode07_initial_sorted.bam
    ├── barcode07_masked_sorted.bam
    ├── barcode07_baseline_sorted.bam
    ├── barcode07.log
    └── ...
```

---

## Usage

### ONT
```bash
./ONT_MITO_v32.sh <input_fastq_folder>
```

### Illumina
```bash
./ILM_MITO_v32.sh <input_fastq_folder>
```

The input must be a directory containing the expected FASTQ files.

- ONT: `.fastq` or `.fastq.gz`
- Illumina: paired-end files containing `R1/R2` or `.R1/.R2`

---

## Required Files

### 1. Adapter file
Required for cutadapt trimming.

Default expected filename:

```bash
Updated_Adapter_Primer_List_Cutadapt_cleaned.txt
```

Manual override:

```bash
export ADAPTER_FILE=/path/to/Updated_Adapter_Primer_List_Cutadapt_cleaned.txt
```

### 2. Reference FASTA
Required linearized mitochondrial reference:

```bash
export ref="linearized_mtdna.fasta"
```

### 3. Regions BED
Required annotation BED file:

```bash
export regions_bed="linearized_regions.bed"
```

### 4. ONT minimap2 index
For ONT, the minimap2 index is expected as:

```bash
export mmi_index="linearized_mtdna.fasta.mmi"
```

If missing, it is generated automatically.

---

## Configurable Parameters

The pipelines can be customized by exporting environment variables before execution.

### Common parameters

| Variable | Default | Description |
| --- | --- | --- |
| `threads` | 8 | Number of threads |
| `MIN_LEN` | 90 | Minimum read length after trimming |
| `EXTRA_TRIM` | 0 | Fixed number of bases trimmed from both ends |
| `CUTADAPT_ERR` | 0.10 | Cutadapt error rate |
| `CUTADAPT_OVL` | 5 | Minimum cutadapt overlap |
| `ref` | `linearized_mtdna.fasta` | Reference FASTA |
| `regions_bed` | `linearized_regions.bed` | BED used for annotation |
| `SAFETY_QUAL` | 20 | Initial VCF QUAL threshold |
| `STRICT_QUAL` | 60 | Stricter SNP-only QUAL threshold |
| `MIN_DEPTH` | 10 | Minimum depth after calling |
| `PILEUP_MAX_DEPTH` | 100000 | Maximum depth used by mpileup |
| `BASEQ_MIN` | 20 | Minimum base quality for mpileup |
| `MAPQ_MIN` | 20 | Minimum mapping quality for mpileup |
| `EDGE_PROTECT` | 1 | Boundary rescue window for masking |

### ONT-specific parameters

| Variable | Default | Description |
| --- | --- | --- |
| `MAX_LEN` | 1500 | Maximum read length retained before analysis |
| `QS_MIN` | 10 | Minimum QS score from `qs:f:` tag |
| `mmi_index` | `${ref}.mmi` | Minimap2 index path |

### Illumina-specific parameters

| Variable | Default | Description |
| --- | --- | --- |
| `MAX_LEN` | 3000 | Maximum read length retained before merging |
| `READQ` | 20 | fastp qualified quality threshold |
| `UNQUAL_PCT` | 40 | fastp unqualified percent limit |
| `N_BASE_LIMIT` | 5 | Maximum allowed N bases |

---

## Example Usage

### ONT with default settings
```bash
./ONT_MITO_v32.sh /path/to/ont_fastqs
```

### ONT with custom masking boundary rescue
```bash
export EDGE_PROTECT=2
./ONT_MITO_v32.sh /path/to/ont_fastqs
```

### Illumina with custom merge and trimming settings
```bash
export READQ=25
export UNQUAL_PCT=30
export EXTRA_TRIM=5
./ILM_MITO_v32.sh /path/to/illumina_fastqs
```

---

## Software Requirements

The following tools must be available in `PATH`.

### Common
- bash
- fastqc
- cutadapt
- samtools
- bcftools
- python3
- awk

### ONT pipeline
- fastplong
- minimap2

### Illumina pipeline
- fastp
- bwa-mem2

---

## Interpretation Notes

These pipelines are diagnostic by design.

The **masked track** is the intended primary callset, while the **baseline track** and additional masked outputs are retained so that masking behavior can be checked against the unmasked data. This makes it possible to:

- confirm that expected variants are still retained
- identify potential masking-related dropout
- inspect homopolymer and blacklist behavior
- compare depth loss between masked and unmasked tracks
- troubleshoot unexpected changes in final calls

For routine downstream review, the file of primary interest is:

- **`*_masked_snps.vcf`**

---

## Summary

Version 32 is a **dual-track diagnostic mtDNA pipeline** for ONT and Illumina data that produces a masked and unmasked comparison framework around the same input sample.

Its main purpose is to support confident review of a masking-based calling strategy while preserving enough diagnostic output to validate the behavior of that strategy on real data.

## Important Note
Results from the pipeline are based on linearized coordinates. To reveal the true positions of variants, they must be reverted back to their original locations using the post processing files provided.
