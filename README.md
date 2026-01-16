# ONT DualCall SNP and INDEL Pipeline

A reproducible Oxford Nanopore Technologies (ONT) pipeline for mitochondrial DNA (mtDNA) variant calling that generates **two VCF outputs per sample from a single callset**: one SNP-only and one SNP plus INDELs.

---

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Pipeline Logic](#pipeline-logic)
- [Requirements](#requirements)
- [Installation](#installation)
- [Inputs](#inputs)
- [Usage](#usage)
- [Outputs](#outputs)
- [Configuration](#configuration)
- [Repository Structure](#repository-structure)
- [Intended Use](#intended-use)
- [Limitations](#limitations)
- [Reproducibility](#reproducibility)
- [License](#license)

---

## Overview

This repository provides a single Bash-based pipeline for conservative mtDNA variant calling from ONT sequencing data.  
The pipeline performs alignment and variant calling **once per sample** and derives two downstream VCFs from the same filtered callset:

- **SNP-only VCF**
- **SNP + INDEL VCF**

This design guarantees direct comparability between variant representations without introducing caller-level divergence.

---

## Features

- Single-pass variant calling per sample
- Dual-output VCF generation from a shared callset
- Conservative consensus calling strategy
- Transparent preprocessing and filtering steps
- Per-sample logging and QC outputs
- Minimal dependencies, no custom binaries

---

## Pipeline Logic

FASTQ
↓
FastQC
↓
cutadapt (3′ pass)
↓
cutadapt (5′ pass)
↓
Fixed-end trimming
↓
minimap2 (map-ont)
↓
Sorted BAM
↓
bcftools mpileup
↓
bcftools call (consensus, ploidy 1)
↓
QUAL filtering
↓
┌───────────────────────┬────────────────────────┐
│ SNP-only VCF          │ SNP + INDEL VCF         │
│ (INDELs removed)      │ (all variants retained)│
└───────────────────────┴────────────────────────┘

---

## Requirements

All tools must be available on the system PATH.

### Software

| Tool | Purpose |
|-----|---------|
| bash | Pipeline execution |
| fastqc | Read and BAM quality control |
| cutadapt | Adapter and primer trimming |
| minimap2 | ONT read alignment |
| samtools | BAM processing |
| bcftools | Variant calling |

---

## Installation

### Conda (Recommended)

```bash
conda create -n ont_dualcall \
  fastqc cutadapt minimap2 samtools bcftools \
  -c bioconda -c conda-forge

conda activate ont_dualcall


⸻

Inputs

FASTQ Files
	•	Input is a directory containing .fastq files
	•	One FASTQ file per sample

input_fastqs/
├── sample_01.fastq
├── sample_02.fastq


⸻

Reference Genome
	•	Mitochondrial reference genome in FASTA format
	•	Linearized mtDNA reference recommended

export ref=linearized_mtdna.fasta


⸻

Adapter File

Required cutadapt adapter list:

Updated_Adapter_Primer_List_Cutadapt_cleaned.txt

Search order:
	1.	Parent directory of the FASTQ folder
	2.	Current working directory

Optional explicit specification:

export ADAPTER_FILE=/path/to/Updated_Adapter_Primer_List_Cutadapt_cleaned.txt


⸻

Usage

./ont_dualcall_snps_indels_v1.0.sh input_fastqs/


⸻

Outputs

Per sample, the pipeline produces:

File	Description
*_all_variants_raw.vcf	Unfiltered SNP + INDEL calls
*_all_variants.vcf	Filtered SNP + INDEL VCF
*_snps.vcf	SNP-only VCF
*_sorted.bam	Sorted alignment
*_sorted.bam.bai	BAM index
*.log	Execution log
FastQC reports	QC summaries


⸻

Configuration

Parameters may be overridden via environment variables:

export threads=8
export QUAL_MIN=20
export BASEQ_MIN=20
export MAPQ_MIN=20
export MIN_LEN=90
export EXTRA_TRIM=22
export PILEUP_MAX_DEPTH=100000


⸻

Repository Structure

.
├── ont_dualcall_snps_indels_v1.0.sh
├── README.md
├── LICENSE
└── example/
    ├── input_fastqs/
    └── expected_output/


⸻

Intended Use
	•	mtDNA variant analysis
	•	SNP-only vs SNP+INDEL comparison
	•	Method benchmarking and validation
	•	Forensic and research workflows

⸻

Limitations
	•	Consensus calling only
	•	No heteroplasmy modeling
	•	No indel normalization or realignment
	•	Not validated for clinical diagnostics

⸻

Reproducibility
	•	Single callset per sample
	•	Deterministic SNP derivation
	•	Per-sample logs capture all commands
	•	External version pinning recommended


