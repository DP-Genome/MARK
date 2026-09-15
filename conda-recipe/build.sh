#!/bin/bash

# Create the bin directory in the conda installation prefix
mkdir -p $PREFIX/bin

# Copy the Bash pipelines
cp MARK.sh $PREFIX/bin/ || echo "Warning: MARK.sh not found"
cp MARK-I.sh $PREFIX/bin/ || echo "Warning: MARK-I.sh not found"

# Copy the Python GUI and scripts
cp MARKLaunch.py $PREFIX/bin/ || echo "Warning: MARKLaunch.py not found"
cp SequenceAnalyzerApp_v3.py $PREFIX/bin/ || true
cp vcf_compare_guiv3.py $PREFIX/bin/ || true
cp VCF_OrganizerApp.py $PREFIX/bin/ || true

# Copy necessary reference files so the GUI can find them automatically
cp linearized_mtdna.fasta $PREFIX/bin/ || echo "Warning: Fasta missing"
cp linearized_regions.bed $PREFIX/bin/ || echo "Warning: BED missing"
cp MARK_Adapter_List_Illumina.txt $PREFIX/bin/ || echo "Warning: Illumina adapter list missing"
cp MARK_Adapter_List_ONT.txt $PREFIX/bin/ || echo "Warning: ONT adapter list missing"
cp rCRS.fasta $PREFIX/bin/ || echo "Warning: rCRS missing"

# Measured CRM Nested coordinates: amplicon inserts, full PCR products, and the
# primer footprints between them. Reference/annotation data - the pipelines carry
# their own copy of the product coordinates internally.
cp CRM_Nested_inserts.bed $PREFIX/bin/ || echo "Warning: inserts BED missing"
cp CRM_Nested_products.bed $PREFIX/bin/ || echo "Warning: products BED missing"
cp CRM_Nested_primers_empirical.bed $PREFIX/bin/ || echo "Warning: primers BED missing"

# Copy the demo test file
cp Test_M.fastq $PREFIX/bin/ || echo "Warning: Test file missing"

# Ensure the executable scripts have the correct permissions
chmod +x $PREFIX/bin/*.sh
chmod +x $PREFIX/bin/*.py
