# Mitochondrial DNA Pipeline Helper Tools

## Overview

This repository contains four companion programs that support the mitochondrial DNA pipeline workflow beyond the core ONT and Illumina shell pipelines. Together, these tools help with:

- launching and controlling the pipeline from a GUI
- applying post-processing to pipeline outputs
- organizing corrected VCF files into biologically meaningful groups
- comparing VCF concordance across replicates or experiments
- inspecting raw, preprocessed, and aligned sequence files to understand read retention and quality progression

These programs are not standalone replacements for the main pipelines. They are support utilities designed to make pipeline execution, QC, organization, and downstream comparison easier and more transparent.

## Included Programs

- `MitoPipelineDashboardv4.7.py`
- `VCF_OrganizerApp.py`
- `vcf_compare_guiv3.py`
- `SequenceAnalyzerApp_v3.py`

They were documented based directly on the uploaded source files. The dashboard logic comes from `MitoPipelineDashboardv4.7.py`, the VCF comparison logic from `vcf_compare_guiv3.py`, the organization utility from `VCF_OrganizerApp.py`, and the sequence progression analyzer from `SequenceAnalyzerApp_v3.py`.

---

## Suggested Workflow

The tools fit together in the following order:

1. **Run the main pipeline** with `MitoPipelineDashboardv4.7.py`.
2. **Organize raw outputs, correct VCF positions, clean BAM headers, and gather final results** using the dashboard’s post-processing tab or full auto cycle.
3. **Organize corrected VCFs by biological sample or barcode group** using `VCF_OrganizerApp.py`.
4. **Compare concordance across VCFs** using `vcf_compare_guiv3.py`.
5. **Inspect read progression and retention from raw FASTQ to processed FASTQ to BAM** using `SequenceAnalyzerApp_v3.py`.

That is the main logic behind the toolkit. One program runs and standardizes the pipeline outputs, one organizes final VCF products, one compares biological concordance, and one helps explain what happened to the reads during processing.

---

## 1. MitoPipelineDashboardv4.7.py

### Purpose

This is the main control panel for running your shell pipelines and then post-processing their outputs from one interface. It is a Tkinter GUI with two tabs:

- **Step 1: Run Pipeline**
- **Step 2: Post-Processing**

It is the operational hub of the workflow.

### Main logic

The dashboard does three important things:

First, it allows you to run a selected pipeline script while exposing key environment-variable parameters in a user-editable form. It reads defaults from the selected shell script and lets the user override only the values they want.

Second, it tracks run progress, captures the pipeline output folder from the line `All outputs written to: ...`, and supports stopping the entire process by terminating the spawned process group, not just the parent shell. That matters because your pipelines spawn child commands such as `samtools`, `minimap2`, or `bwa-mem2`.

Third, it performs the downstream cleanup needed to make the outputs easier to interpret and compare. That includes collecting VCFs and BAMs, correcting VCF positions from the linearized reference back to original coordinates, cleaning BAM headers, and assembling a final results folder.

### What it does

#### Pipeline tab

The pipeline tab lets you:

- select the shell script to run
- inspect script defaults and optionally override them
- choose the linearized reference FASTA, BED file, and adapter file
- choose the input FASTQ folder
- optionally choose an output base directory
- optionally assign a custom run name
- run the pipeline only
- run the pipeline plus all post-processing automatically
- stop a running job cleanly

The parameter grid includes values such as `MAX_LEN`, `QS_MIN`, `STRICT_QUAL`, `MIN_DEPTH`, `MIN_LEN`, `MIN_LEN_POST`, `MAX_LEN_POST`, `EXTRA_TRIM`, `READQ`, and `UNQUAL_PCT`. The dashboard only passes overrides for variables whose checkbox is enabled.

#### Post-processing tab

The post-processing tab is broken into four steps:

1. **Organize Raw Output**
   - recursively copies VCFs into a `vcfs/` folder and sorted BAMs into a `sorted_bams/` folder inside the run directory.

2. **Linear VCF Correction**
   - converts linearized variant positions back to original circular mtDNA coordinates using a position map built from the original reference. Output files are named with `_corrected.vcf`.

3. **BAM Header Cleaning**
   - uses `samtools view -H`, `sed`, and `samtools reheader` to produce `_cleaned.bam` files and new BAM indexes. This is intended to improve downstream compatibility, especially for viewing or reusing the BAMs.

4. **Final Collection**
   - gathers corrected VCFs and cleaned BAMs into a single `Final_Pipeline_Results/` directory with `corrected_vcfs/` and `cleaned_bams/` subfolders.

### How to use it

Run the GUI:

```bash
python3 MitoPipelineDashboardv4.7.py
```

Typical use:

1. Open the program.
2. In the first tab, browse to the pipeline shell script you want to run.
3. Set the reference FASTA, BED file, adapter list, and input FASTQ folder.
4. Leave parameters at script defaults unless you want to override specific values.
5. Click **RUN PIPELINE ONLY** if you want only the shell pipeline.
6. Click **RUN PIPELINE & POST-PROCESSING (AUTO)** if you want the full workflow in one pass.

### Important notes

- The default script path inside the program points to `ONT_MITO_v36.sh`, so for newer pipeline versions you should use **Browse & Load** and select the correct `.sh` file manually.
- The auto cycle is the most practical mode if you want one-button execution from raw pipeline run to cleaned final outputs.
- The VCF coordinate correction assumes the mtDNA split position is `8284`.

---

## 2. VCF_OrganizerApp.py

### Purpose

This program organizes final corrected VCF files into grouped folders so they are easier to review by biological sample, experiment, or ad hoc search terms. It is a file-management tool, not a variant-analysis engine.

### Main logic

The organizer searches recursively for VCF files matching a chosen suffix, then copies them into a new collection folder. It supports two modes:

- **Auto Organization** using a sample mapping JSON
- **Manual Extraction / Custom Search** using filename substring matching

The default target suffix is:

```text
_trimmed_snps_corrected.vcf
```

That default reflects the current trimmed-track primary output pattern after correction.

### What it does

#### Auto Organization

This mode loads a JSON mapping file in which each biological sample name maps to one or more barcode or filename identifiers. It then scans all matching VCFs and copies each file into the correct sample folder if the filename contains any identifier associated with that sample.

Expected mapping structure example:

```json
{
  "2800M": ["barcode01", "barcode02", "barcode03"],
  "007": ["barcode04", "barcode05", "barcode06"]
}
```

The exact strings can be barcodes, sample tokens, or any other consistent filename markers.

#### Manual Extraction

This mode is simpler. You provide comma-separated search substrings and a target folder name, and the program copies every matching VCF into that folder. This is useful for extracting a subset of files quickly without changing the main mapping.

### How to use it

Run the GUI:

```bash
python3 VCF_OrganizerApp.py
```

#### Auto organization workflow

1. Set the input directory containing pipeline outputs or final corrected VCFs.
2. Set the destination output directory.
3. Select the sample mapping JSON.
4. Keep or change the VCF suffix filter.
5. Optionally enter a master folder base name.
6. Click **Run Auto Organize**.

#### Manual extraction workflow

1. Set input and output directories.
2. Enter search substrings separated by commas.
3. Enter the name of the folder to create.
4. Click **Search & Extract**.

### Output behavior

The tool creates a timestamped master collection folder and then creates per-sample or per-query subfolders inside it. It copies files rather than moving them, which is the correct choice for preserving the original pipeline results.

---

## 3. vcf_compare_guiv3.py

### Purpose

This is the concordance and VCF comparison tool. It reads multiple VCF files, extracts variants, determines which sites are shared or discordant, and generates both tabular and visual summaries. It is intended for comparing replicates, comparing pipeline versions, or checking agreement between callsets.

### Main logic

The core logic is based on a variant presence matrix.

For each VCF, the program parses the variant records, extracts basic fields such as chromosome, position, REF, ALT, QUAL, DP, type, and inferred allele frequency where available, then builds a combined table keyed by variant position.

Each variant site is then classified as one of three categories:

- **Consensus (All)**: present in all selected VCFs
- **Discordant (Some)**: present in more than one file but not all
- **Unique (Singleton)**: present in only one file

The program also annotates whether a position falls inside predefined homopolymer regions and whether the variant is a SNP or INDEL.

### What it produces

For each run, it creates a timestamped results folder containing:

- `variants_summary_matrix.csv`
- `variants_detailed_raw.csv`
- `file_statistics.csv`
- `plot_presence_by_position.png`
- `plot_problematic_variants.png`
- `plot_quality_heatmap.png`

### Plot logic

#### Presence plot

Shows whether each variant site is present or absent across files and overlays depth values inside the colored cells. Regions HV1, HV2, HV3, and Coding are visually marked.

#### Problematic variants plot

Highlights INDELs and homopolymer-associated SNPs, separating them visually from cleaner SNPs. This is useful when checking whether discordance clusters in known difficult regions.

#### Quality heatmap

Colors each present variant call by quality score and writes the actual quality value in the cell. This makes it easier to distinguish a true disagreement from a weak-call disagreement.

### How to use it

Run the GUI:

```bash
python3 vcf_compare_guiv3.py
```

Typical workflow:

1. Add individual VCF files or a folder containing VCFs.
2. Choose an output directory.
3. Select at least two VCFs.
4. Click **RUN ANALYSIS**.

### When to use it

Use this tool when you want to:

- compare triplicates
- compare ONT and Illumina outputs
- compare baseline vs trimmed pipeline outputs
- compare two pipeline versions
- find which variants are unique, shared, or missing across replicates
- identify whether disagreements are enriched in homopolymer regions or low-quality regions

---

## 4. SequenceAnalyzerApp_v3.py

### Purpose

This program is a bulk sequence-file analyzer for BAM and FASTQ files. Its job is to help explain what happened to the reads across processing stages by comparing raw, preprocessed, and aligned sequence metrics.

It is a QC and progression-inspection tool rather than a variant caller.

### Main logic

The program scans a pipeline run folder, identifies sample directories, and looks for representative files from three stages:

- a raw-stage FASTQ
- a preprocessed FASTQ
- an aligned BAM

It then calculates metrics such as:

- read count
- average read length
- read-length distribution
- average base quality
- mapping percentage for BAM files

The idea is to show how the data changed between stages, not just whether the final VCF exists.

### What it produces

For each sample, it generates progression plots showing:

- **Read Length Progression**
- **Base Quality Progression**

At the run level it also generates:

- a TSV summarizing length-distribution bins across stages
- an optional offline HTML dashboard that lets the user inspect all samples in one place.

### Sample-stage detection logic

The program tries to infer representative files using filename patterns. For example, it looks for raw-stage files such as `*_lenLE1500.fastq`, `*_lenLE3000_R1.fastq`, `*_qsGE10.fastq`, or `*_merged.fastq`, then preprocessing-stage files such as `*_trim5.fastq` or `*_qsGE10_trim5.fastq`, and finally aligned BAMs such as `*trimmed_sorted.bam` or `*initial_sorted.bam`.

That means the tool is tightly coupled to your mitochondrial pipeline naming conventions. This is a strength when the pipeline outputs follow those conventions, and a limitation when they do not.

### How to use it

Run the GUI:

```bash
python3 SequenceAnalyzerApp_v3.py
```

Typical workflow:

1. Select the pipeline output folder containing per-sample subdirectories.
2. Select an output folder.
3. Enter a run name or click **Auto**.
4. Leave HTML dashboard generation enabled unless you specifically do not want it.
5. Click **Run Analysis**.

### When to use it

Use this tool when you want to understand:

- where reads are being lost
- whether trimming shifted the read-length profile too aggressively
- whether BAM read quality or mapping fraction looks abnormal
- whether a problematic sample looks obviously different from the others before you start analyzing concordance

---

## End-to-End Practical Use Case

A practical use case for the full set of tools would look like this:

1. Run `ONT_MITO_v36.sh` or `ILM_MITO_v36.sh` through `MitoPipelineDashboardv4.7.py`.
2. Let the dashboard auto-run post-processing so you end up with corrected VCFs and cleaned BAMs in `Final_Pipeline_Results/`.
3. Use `VCF_OrganizerApp.py` to group corrected trimmed-track VCFs by biological sample such as `2800M`, `007`, or other sample identities.
4. Use `vcf_compare_guiv3.py` to compare triplicates or cross-platform callsets and identify consensus, discordant, and singleton variants.
5. If a sample behaves abnormally, run `SequenceAnalyzerApp_v3.py` on the pipeline output folder to see whether the issue already appears in the read-length, quality, or mapping progression.

That combined workflow covers execution, cleanup, organization, concordance, and diagnostic interpretation.

---

## Dependencies

These programs rely on a mixture of Python libraries and external command-line tools.

### Python libraries

Depending on which tool you run, you may need:

- `tkinter`
- `Biopython`
- `numpy`
- `pandas`
- `matplotlib`
- `pysam`

Install them with:

```bash
pip install biopython numpy pandas matplotlib pysam
```

### External command-line tools

The dashboard’s post-processing relies on tools such as:

- `samtools`
- `sed`
- the shell pipeline itself and its dependencies

---

## Known Limitations and Notes

- `MitoPipelineDashboardv4.7.py`  defaults to the current script released path that are assumed to be in the same directory as the program itself. If using a different script, it should be pointed to manually.
- `VCF_OrganizerApp.py` depends entirely on filename conventions and a correct mapping JSON. If the sample identifiers in filenames are inconsistent, organization will be unreliable.
- `vcf_compare_guiv3.py` compares by parsed VCF position and alleles. It is useful for concordance analysis, but it is not a truth-evaluation engine with sensitivity and specificity metrics.

---

## Summary

These four helper programs form the supporting software layer around the mitochondrial pipeline.

- `MitoPipelineDashboardv4.7.py` runs the pipeline and standardizes outputs.
- `VCF_OrganizerApp.py` groups corrected VCFs into meaningful collections.
- `vcf_compare_guiv3.py` measures agreement and disagreement across VCF files.
- `SequenceAnalyzerApp_v3.py` explains read-level progression from raw data to alignment.

Taken together, they give the user a practical way to operate the pipeline, clean the outputs, organize the result files, inspect biological concordance, and troubleshoot sample behavior.
