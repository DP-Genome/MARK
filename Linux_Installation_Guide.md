# MARK Pipeline - Linux Installation and Usage Guide

This guide will walk you through setting up and running the **MARK (Mitochondrial Amplicon Resolving Kit)** pipeline on a Linux system.

## 1. Installation

The easiest way to install the pipeline and all its dependencies is using Conda. We recommend creating a fresh environment to avoid conflicts.

Open your terminal and run the following commands:

```bash
# 1. Create a new conda environment with the mark_pipeline package, in one command
#    (creating the environment first and installing into it can make conda pick an old MARK version)
conda create -n mark_env -c akmartian -c conda-forge -c bioconda mark_pipeline -y

# 2. Install freebayes, MARK's variant caller, in its own environment named mark_fb (required).
#    Its builds need an older htslib than MARK's other tools, so they cannot share one environment.
#    MARK finds the mark_fb environment by itself; nothing to configure.
conda create -n mark_fb -c conda-forge -c bioconda freebayes=1.3.10 -y

# 3. Activate MARK's environment
conda activate mark_env

# 4. Check the version (MARK.sh --help prints it on the first line)
MARK.sh --help | head -1
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

**Note:** The core pipelines (`MARK.sh` and `MARK-I.sh`) are fully standalone scripts. They can be executed natively in the terminal without a graphical interface. This is ideal for HPC environments, workflow automation, or overriding parameters manually. If you prefer a visual interface, you can optionally launch the GUI dashboard using `MARKLaunch.py`.

Here is a quick example of how to run the pipeline on your FastQ data using the command line.

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
