#!/bin/bash
# collect_corrected_vcfs_and_cleaned_bams.sh
# Usage: ./collect_corrected_vcfs_and_cleaned_bams.sh /path/to/main_folder

set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <main_folder>"
  exit 1
fi

main_dir="$1"

if [ ! -d "$main_dir" ]; then
  echo "Error: '$main_dir' is not a directory."
  exit 1
fi

# Source folders (from the first script)
vcf_src_dir="$main_dir/vcfs"
bam_src_dir="$main_dir/sorted_bams"

if [ ! -d "$vcf_src_dir" ]; then
  echo "Warning: VCF source directory '$vcf_src_dir' not found. Skipping VCF search."
fi

if [ ! -d "$bam_src_dir" ]; then
  echo "Warning: BAM source directory '$bam_src_dir' not found. Skipping BAM search."
fi

# Destination folders
clean_vcf_dir="$main_dir/cleaned vcfs"
clean_bam_dir="$main_dir/cleaned bams"

mkdir -p "$clean_vcf_dir" "$clean_bam_dir"

echo "Main directory:       $main_dir"
echo "VCF source dir:       $vcf_src_dir"
echo "BAM source dir:       $bam_src_dir"
echo "Corrected VCFs dir:   $clean_vcf_dir"
echo "Cleaned BAMs dir:     $clean_bam_dir"
echo

########################################
# Helper: copy file with collision-safe name
########################################
copy_with_suffix() {
  src="$1"
  dest_dir="$2"

  base="$(basename "$src")"
  name="${base%.*}"
  ext="${base##*.}"

  # Handle files without an extension
  if [ "$name" = "$base" ]; then
    name="$base"
    ext=""
  fi

  if [ -n "$ext" ]; then
    candidate="$dest_dir/$base"
  else
    candidate="$dest_dir/$name"
  fi

  if [ -e "$candidate" ]; then
    i=1
    if [ -n "$ext" ]; then
      while [ -e "$dest_dir/${name}_$i.$ext" ]; do
        i=$((i+1))
      done
      candidate="$dest_dir/${name}_$i.$ext"
    else
      while [ -e "$dest_dir/${name}_$i" ]; do
        i=$((i+1))
      done
      candidate="$dest_dir/${name}_$i"
    fi
  fi

  echo "Copying: $src"
  echo "   to  : $candidate"
  cp "$src" "$candidate"
}

########################################
# 1. Collect corrected VCF files
########################################
if [ -d "$vcf_src_dir" ]; then
  echo "Searching for corrected VCF files in '$vcf_src_dir'..."

  find "$vcf_src_dir" -type f \( \
      -iname '*corrected*.vcf'    -o \
      -iname '*corrected*.vcf.gz' \
    \) -print0 |
  while IFS= read -r -d '' file; do
    copy_with_suffix "$file" "$clean_vcf_dir"
  done

  echo "Finished collecting corrected VCF files."
  echo
fi

########################################
# 2. Collect cleaned BAMs and index files (still using 'cleaned' in name)
########################################
if [ -d "$bam_src_dir" ]; then
  echo "Searching for cleaned BAMs and index files in '$bam_src_dir'..."

  find "$bam_src_dir" -type f \( \
      -iname '*cleaned*.bam'      -o \
      -iname '*cleaned*.bam.bai'  -o \
      -iname '*cleaned*.bai' \
    \) -print0 |
  while IFS= read -r -d '' file; do
    copy_with_suffix "$file" "$clean_bam_dir"
  done

  echo "Finished collecting cleaned BAMs and indexes."
  echo
fi

echo "Done."
