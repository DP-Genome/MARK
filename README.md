# ONT_MITO_CALL_version2

## Overview

This repository contains a specialized bioinformatics pipeline for the analysis of mitochondrial DNA (mtDNA) sequencing data generated using Oxford Nanopore Technologies (ONT).

The pipeline accepts **either a single FASTQ file or multiple FASTQ files contained within a directory**. Each FASTQ file is processed independently.

The pipeline performs **read quality filtering (QS score), rigorous trimming, alignment, variant calling, and region-based annotation**. Unlike standard pipelines, this workflow produces **four distinct VCF datasets per sample** to separate high-confidence variants from homopolymer-associated calls and known artifacts.

Correct operation depends on:

* A **supplied adapter/primer/barcode file** for trimming.
* A **linearized mitochondrial DNA reference genome**.
* A **BED file** defining Homopolymer (HP) and Blacklisted regions for annotation.

Core tools include `cutadapt`, `minimap2`, `samtools`, and `bcftools`.

---

## Pipeline Workflow

The major processing steps are described below.

---

## Step 1: Quality Score (QS) Filtering

Before processing, raw reads are filtered based on the `qs:f` tag found in standard ONT FASTQ headers. This ensures only high-quality reads enter the analysis pipeline.

* **Default Threshold:** Q-Score ≥ 10
* *Configurable via `QS_MIN` environment variable.*

```bash
# Internal logic
awk '... keep read if qs >= 10 ...' input.fastq > input_qsGE10.fastq

```

---

## Step 2: Quality Control (FastQC)

Quality assessment is performed on the **QS-filtered** reads to verify improvement and identify remaining issues.

```bash
fastqc sample_qsGE10.fastq

```

---

## Step 3: Adapter & Primer Trimming (cutadapt)

Adapters and primers are removed using a **two-pass trimming strategy**:

1. **3′ adapter trimming**
2. **5′ adapter trimming**
3. **Fixed End Trimming:** A fixed number of bases (default: 22) are trimmed from *both* ends of the read to remove residual artifacts.
4. **Length Filtering:** Reads shorter than the minimum length (default: 90bp) are discarded.

```bash
cutadapt -u 22 -u -22 -m 90 ...

```

---

## Step 4: Sequence Alignment (minimap2)

Trimmed reads are aligned to a **linearized mitochondrial DNA reference genome** using the `map-ont` preset.

```bash
minimap2 -ax map-ont linearized_mtdna.fasta trimmed_reads.fastq | samtools sort > sorted.bam

```

---

## Step 5: Variant Calling (bcftools)

Variant calling is performed using a consensus-based approach (Ploidy 1).

* **High Depth Support:** `mpileup` is configured to handle depths up to 100,000x.
* **Qual Filtering:** Low-confidence variants (QUAL < 20) are discarded immediately.

```bash
bcftools mpileup -d 100000 -Q 20 -q 20 ... | bcftools call -cv --ploidy 1

```

---

## Step 6: Region Annotation & Output Splitting

This is the core differentiation of the v2 pipeline. The raw variants are **annotated** using a user-provided BED file to flag regions as either "HP_Region" (Homopolymer) or "Blacklist_Site".

The annotated VCF is then split into four specific output files:

1. **Annotated All:** All variants with region tags.
2. **SNPs (Special):** SNPs only. **Excludes** known artifact positions (e.g., 7898, 7899, 8595).
3. **Clean:** High-confidence variants. **Excludes** all Homopolymer regions and Blacklist sites.
4. **Homopolymers:** Contains **only** variants found within Homopolymer regions.

---

## Outputs

For each input FASTQ file (e.g., `SampleA`), the pipeline creates a folder `ONT_MITO_CALL_version2_[TIMESTAMP]_output/SampleA/` containing:

| File Name | Description |
| --- | --- |
| `*_qsGE10_clean.vcf` | **Primary Output.** Variants excluding HP regions and Blacklist sites. |
| `*_qsGE10_snps.vcf` | SNPs only. Excludes specific artifact positions (7898, 7899, 8595). |
| `*_qsGE10_homopolymers.vcf` | Variants found exclusively in Homopolymer regions. |
| `*_qsGE10_annotated_all.vcf` | The complete callset with `RegionType` annotations. |
| `*_qsGE10.fastq` | The QS-filtered FASTQ file used for analysis. |
| `*.bam` / `*.bam.bai` | Sorted alignment file and index. |
| `*.log` | Detailed execution log. |

---

## Usage

```bash
./ONT_MITO_CALL_version2.sh <input_fastq_folder>

```

The input must be a directory containing `.fastq` or `.fastq.gz` files.

---

## Configuration & Requirements

### 1. Adapter File

A cutadapt adapter file is required. The pipeline looks for `Updated_Adapter_Primer_List_Cutadapt_cleaned.txt` in the parent or current directory, or you can specify it manually:

```bash
export ADAPTER_FILE=/path/to/adapters.txt

```

### 2. Reference & Annotation Files

You must provide the linearized reference and the corresponding regions BED file.

```bash
export ref="linearized_mtdna.fasta"
export regions_bed="linearized_regions.bed"

```

*(Note: Ensure `regions_bed` is generated using the associated Python script before running the pipeline).*

### 3. Environment Variables (Optional Overrides)

You can override these defaults by exporting them before running the script:

| Variable | Default | Description |
| --- | --- | --- |
| `QS_MIN` | **10** | Minimum Read Mean Quality Score (from `qs:f` tag). |
| `threads` | 8 | Number of threads for parallel processing. |
| `QUAL_MIN` | 20 | Minimum Variant QUAL score. |
| `MIN_LEN` | 90 | Minimum read length after trimming. |
| `EXTRA_TRIM` | 22 | Bases trimmed from both ends of the read. |
| `PILEUP_MAX_DEPTH` | 100000 | Max depth for pileup (prevents downsampling). |

---

## Software Requirements

The following tools must be in your PATH:

* bash
* fastqc
* cutadapt
* minimap2
* samtools
* bcftools
* awk
