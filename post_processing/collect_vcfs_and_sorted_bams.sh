#!/bin/bash
# collect_vcfs_and_sorted_bams.sh
# Usage: ./collect_vcfs_and_sorted_bams.sh /path/to/main_folder

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

# Output folders
vcf_dir="$main_dir/vcfs"
bam_dir="$main_dir/sorted_bams"

mkdir -p "$vcf_dir" "$bam_dir"

echo "Main directory: $main_dir"
echo "VCF output dir: $vcf_dir"
echo "BAM output dir: $bam_dir"
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

  # For files without an extension (unlikely here), handle gracefully
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
# 1. Collect VCF files
########################################
echo "Searching for VCF files..."

find "$main_dir" \
  \( -path "$vcf_dir" -o -path "$bam_dir" \) -prune -o \
  -type f \( -iname '*.vcf' -o -iname '*.vcf.gz' \) -print0 |
while IFS= read -r -d '' file; do
  copy_with_suffix "$file" "$vcf_dir"
done

echo "Finished collecting VCF files."
echo

########################################
# 2. Collect sorted BAMs and indexes
########################################
echo "Searching for sorted BAMs and index files..."

find "$main_dir" \
  \( -path "$vcf_dir" -o -path "$bam_dir" \) -prune -o \
  -type f \( \
    -iname '*sorted.bam'       -o \
    -iname '*sorted.bam.bai'   -o \
    -iname '*sorted.bai' \
  \) -print0 |
while IFS= read -r -d '' file; do
  copy_with_suffix "$file" "$bam_dir"
done

echo "Finished collecting sorted BAMs and indexes."
echo "Done."
