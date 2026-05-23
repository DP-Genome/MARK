# MARK Pipeline - Linux Installation and Usage Guide

This guide will walk you through setting up and running the **MARK (Mitochondrial Amplicon Resolving Kit)** pipeline on a Linux system.

## 1. Installation

The easiest way to install the pipeline and all its dependencies is using Conda. We recommend creating a fresh environment to avoid conflicts.

Open your terminal and run the following commands:

```bash
# 1. Create a new conda environment for the pipeline
conda create -n mark_env python=3.10 -y

# 2. Activate the new environment
conda activate mark_env

# 3. Install the mark_pipeline package
conda install -c akmartian -c conda-forge -c bioconda mark_pipeline -y
```

## 2. Verify Installation

Once installed, verify that the pipeline tools are properly loaded in your path by running:

```bash
MARK.sh --help
MARK-I.sh --help
```
*(Both commands should print their usage instructions without errors).*

If you are using a system with a graphical interface (GUI), you can also verify the launcher:
```bash
MARKLaunch.py
```

## 3. Usage Example (ONT Data Test)

Here is a quick example of how to run the pipeline on your FastQ data.

```bash
# 1. Create a working directory and an input folder for your FASTQ files
mkdir -p ~/mark_analysis/input_fastqs
cd ~/mark_analysis

# 2. Copy or move your sample FASTQ files into the input folder
cp /path/to/your/sample.fastq.gz ~/mark_analysis/input_fastqs/

# 3. Run the pipeline and save the output log
MARK.sh ~/mark_analysis/input_fastqs 2>&1 | tee MARK_run.log
```

### Checking Your Results
You can monitor the tail of your log file to see the progress:
```bash
tail -n 100 MARK_run.log
```
Once completed, all output directories and files will be generated inside the `~/mark_analysis` folder.
