#!/bin/bash

# Create the bin directory in the conda installation prefix
mkdir -p $PREFIX/bin

# Copy the Bash pipelines
cp MARK.sh $PREFIX/bin/ || echo "Warning: MARK.sh not found"
cp MARK-I.sh $PREFIX/bin/ || echo "Warning: MARK-I.sh not found"

# Copy the Python GUI and scripts
cp MARKLaunch.py $PREFIX/bin/ || echo "Warning: MARKLaunch.py not found"
cp post_processing/SequenceAnalyzerApp_v3.py $PREFIX/bin/ || true
cp post_processing/vcf_compare_guiv3.py $PREFIX/bin/ || true
cp post_processing/VCF_OrganizerApp.py $PREFIX/bin/ || true

# Copy necessary reference files so the GUI can find them automatically
cp linearized_mtdna.fasta $PREFIX/bin/ || echo "Warning: Fasta missing"
cp linearized_regions.bed $PREFIX/bin/ || echo "Warning: BED missing"
cp Updated_Adapter_Primer_List_Cutadapt_cleaned.txt $PREFIX/bin/ || echo "Warning: Adapter missing"
cp rCRS.fasta $PREFIX/bin/ || echo "Warning: rCRS missing"

# Copy the demo test file
cp Test_M.fastq $PREFIX/bin/ || echo "Warning: Test file missing"

# Ensure the executable scripts have the correct permissions
chmod +x $PREFIX/bin/*.sh
chmod +x $PREFIX/bin/*.py
