ONT DualCall SNP and INDEL Pipeline

Overview

This repository contains a single Oxford Nanopore Technologies (ONT) pipeline for mitochondrial DNA (mtDNA) variant calling that produces two VCF outputs per sample from one unified callset:
	1.	SNP-only VCF
	2.	SNP + INDEL VCF

Both outputs are derived from the same alignment and variant calling step, ensuring direct comparability between SNP-only and SNP+INDEL results.

The pipeline is intentionally conservative, transparent, and reproducible. It is designed for mtDNA analysis workflows where clarity, consistency, and auditability are priorities.

⸻

Pipeline Name

Script

ont_dualcall_snps_indels_v1.0.sh

Tag

ONT_DualCall_SNPs_INDELs_v1p0


⸻

What This Pipeline Does

For each FASTQ file, the pipeline performs the following steps:

1. Quality Control
	•	FastQC on raw FASTQ input
	•	FastQC on final aligned BAM

2. Adapter and Primer Trimming
	•	Two-pass cutadapt trimming using a supplied adapter list:
	•	3′ adapter trimming
	•	5′ adapter trimming
	•	Additional fixed trimming from both ends after adapter removal
	•	Minimum read length enforced after trimming

3. Alignment
	•	Alignment to a mitochondrial reference genome using minimap2 (ONT preset)
	•	Sorted and indexed BAM output

4. Variant Calling
	•	Variant calling performed once per sample using:
	•	bcftools mpileup with base and mapping quality thresholds
	•	bcftools call (consensus caller, ploidy 1)
	•	QUAL-based filtering

5. Dual VCF Generation

From the same filtered callset, the pipeline produces:
	•	SNP-only VCF
INDELs are removed after variant calling.
	•	SNP + INDEL VCF
All variant types are retained.

⸻

Outputs

For each sample (one FASTQ file), the following files are produced:

File	Description
*_all_variants_raw.vcf	Raw SNP + INDEL calls before QUAL filtering
*_all_variants.vcf	Filtered SNP + INDEL VCF
*_snps.vcf	SNP-only VCF derived from the filtered callset
*_sorted.bam	Sorted alignment
*_sorted.bam.bai	BAM index
*.log	Per-sample execution log
FastQC reports	QC summaries

Each sample is processed independently in its own output folder.

⸻

Intended Use

This pipeline is suitable for:
	•	Mitochondrial DNA variant analysis
	•	Reproducible SNP reporting
	•	Comparative analysis of SNP-only vs SNP+INDEL callsets
	•	Method benchmarking and validation
	•	Forensic and research workflows requiring conservative calling

This pipeline is not intended for:
	•	Somatic variant detection
	•	Heteroplasmy quantification beyond basic inspection
	•	Advanced indel normalization or realignment strategies
	•	Clinical diagnostics without independent validation

⸻

Requirements

All tools must be installed and available in the system PATH.

Required Software

Tool	Purpose
bash	Pipeline execution
fastqc	Quality control
cutadapt	Adapter and primer trimming
minimap2	ONT read alignment
samtools	BAM processing
bcftools	Variant calling

Recommended Installation (Conda)

conda create -n ont_dualcall \
  fastqc cutadapt minimap2 samtools bcftools \
  -c bioconda -c conda-forge

conda activate ont_dualcall


⸻

Input Requirements

FASTQ Files
	•	Input must be a directory containing FASTQ files
	•	File extension must be .fastq
	•	One sample per FASTQ file

Example:

input_fastqs/
├── sample1.fastq
├── sample2.fastq


⸻

Reference Genome
	•	A mitochondrial reference genome in FASTA format
	•	A linearized mtDNA reference is recommended

Set via environment variable:

export ref=linearized_mtdna.fasta


⸻

Adapter File

The pipeline requires a cutadapt adapter list file named:

Updated_Adapter_Primer_List_Cutadapt_cleaned.txt

The pipeline will search for this file in:
	1.	The parent directory of the input FASTQ folder
	2.	The current working directory

You may also specify it explicitly:

export ADAPTER_FILE=/path/to/Updated_Adapter_Primer_List_Cutadapt_cleaned.txt


⸻

Running the Pipeline

./ont_dualcall_snps_indels_v1.0.sh input_fastqs/


⸻

Output Directory Structure

Each run generates a timestamped output directory:

ont_dualcall_snps_indels_v1p0_YYYYMMDD_HHMMSS_output/
└── sample1/
    ├── sample1_all_variants_raw.vcf
    ├── sample1_all_variants.vcf
    ├── sample1_snps.vcf
    ├── sample1_sorted.bam
    ├── sample1_sorted.bam.bai
    ├── sample1.log
    ├── FastQC reports


⸻

Configuration Options

All parameters can be overridden using environment variables:

export threads=8
export QUAL_MIN=20
export BASEQ_MIN=20
export MAPQ_MIN=20
export MIN_LEN=90
export EXTRA_TRIM=22
export PILEUP_MAX_DEPTH=100000

Default values are conservative and suitable for mtDNA sequencing data.
