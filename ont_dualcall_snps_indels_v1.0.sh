#!/opt/homebrew/bin/bash
# ont_dualcall_snps_indels_v1.0.sh
# tag: ONT_DualCall_SNPs_INDELs_v1p0
#
# PURPOSE
# ---------------------------------------------------------------------------
# One pipeline run produces TWO per-sample VCFs from a single variant callset:
#
# 1) SNP-only VCF
#    - Consensus variant calling
#    - Quality filtering
#    - INDELs removed AFTER calling
#    Output: *_snps.vcf
#
# 2) SNP + INDEL VCF
#    - Same exact callset and filtering
#    - SNPs and INDELs retained
#    Output: *_all_variants.vcf
#
# This pipeline is intentionally conservative and transparent.
# It is designed for reproducible mitochondrial DNA analysis,
# benchmarking, and validation workflows.
#
# PREPROCESSING
# ---------------------------------------------------------------------------
# - FastQC (raw FASTQ)
# - cutadapt 2-pass adapter trimming (3′ then 5′)
# - extra fixed trimming (EXTRA_TRIM) from each end AFTER adapter trimming
# - minimap2 alignment (ONT preset)
# - sorted and indexed BAM
#
# VARIANT CALLING
# ---------------------------------------------------------------------------
# - bcftools mpileup with base and mapping quality thresholds
# - bcftools call (consensus caller, ploidy 1)
# - QUAL-based variant filtering
#
# USAGE
# ---------------------------------------------------------------------------
#   ./ont_dualcall_snps_indels_v1.0.sh <input_fastq_folder>
#
# ENV OVERRIDES
# ---------------------------------------------------------------------------
#   export ref=linearized_mtdna.fasta
#   export threads=8
#   export QUAL_MIN=20
#   export MIN_LEN=90
#   export EXTRA_TRIM=22
#   export PILEUP_MAX_DEPTH=100000
#   export BASEQ_MIN=20
#   export MAPQ_MIN=20
# ---------------------------------------------------------------------------

set -euo pipefail
set -o pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <input_fastq_folder>"
  exit 1
fi

input_folder="$1"
[[ -d "$input_folder" ]] || { echo "Error: Input folder '$input_folder' not found."; exit 1; }

threads="${threads:-8}"
pipeline_name="ont_dualcall_snps_indels_v1p0"

ref="${ref:-linearized_mtdna.fasta}"
mmi_index="${mmi_index:-${ref}.mmi}"

QUAL_MIN="${QUAL_MIN:-20}"
MIN_LEN="${MIN_LEN:-90}"
EXTRA_TRIM="${EXTRA_TRIM:-22}"

PILEUP_MAX_DEPTH="${PILEUP_MAX_DEPTH:-100000}"
BASEQ_MIN="${BASEQ_MIN:-20}"
MAPQ_MIN="${MAPQ_MIN:-20}"

# Adapter file detection
ADAPTER_FILE="${ADAPTER_FILE:-}"
if [[ -z "$ADAPTER_FILE" ]]; then
  main_dir="$(dirname "$input_folder")"
  if [[ -f "$main_dir/Updated_Adapter_Primer_List_Cutadapt_cleaned.txt" ]]; then
    ADAPTER_FILE="$main_dir/Updated_Adapter_Primer_List_Cutadapt_cleaned.txt"
  elif [[ -f "Updated_Adapter_Primer_List_Cutadapt_cleaned.txt" ]]; then
    ADAPTER_FILE="$(pwd)/Updated_Adapter_Primer_List_Cutadapt_cleaned.txt"
  fi
fi
if [[ -z "$ADAPTER_FILE" || ! -f "$ADAPTER_FILE" ]]; then
  echo "Error: adapter file 'Updated_Adapter_Primer_List_Cutadapt_cleaned.txt' not found"
  exit 1
fi

run_ts="$(date +%Y%m%d_%H%M%S)"
run_out="${pipeline_name}_${run_ts}_output"
mkdir -p "$run_out"

run_log() {
  {
    printf 'Running: %s\n' "$1"
    bash -lc "$1"
    printf '\n'
  } 2>&1 | tee -a "$log_file"
}

# Reference indexing
[[ -f "$mmi_index" ]] || minimap2 -d "$mmi_index" "$ref"
[[ -f "${ref}.fai" ]] || samtools faidx "$ref"

shopt -s nullglob
found_any=false

for fq in "$input_folder"/*.fastq; do
  [[ -f "$fq" ]] || continue
  found_any=true

  base="$(basename "$fq" .fastq)"
  sample_out="$run_out/${base}"
  mkdir -p "$sample_out"
  log_file="$sample_out/${base}.log"

  echo "Processing $base" | tee -a "$log_file"
  echo "[Info] Reference: $ref" | tee -a "$log_file"
  echo "[Info] mpileup: -d $PILEUP_MAX_DEPTH -Q$BASEQ_MIN -q$MAPQ_MIN" | tee -a "$log_file"
  echo "[Info] Extra trim: ${EXTRA_TRIM} bp each end (post-adapter trimming)" | tee -a "$log_file"

  # QC raw
  run_log "fastqc '$fq' -o '$sample_out'"

  # cutadapt 2-pass
  t3="$sample_out/${base}_trim3.fastq"
  t5="$sample_out/${base}_trim5.fastq"
  run_log "cutadapt -a file:'$ADAPTER_FILE' --error-rate 0.10 --overlap 5 --minimum-length $MIN_LEN --cores $threads -o '$t3' '$fq' > '$sample_out/${base}_cutadapt_a.log'"
  run_log "cutadapt -g file:'$ADAPTER_FILE' --error-rate 0.10 --overlap 5 --minimum-length $MIN_LEN --cores $threads -o '$t5' '$t3' > '$sample_out/${base}_cutadapt_g.log'"

  # extra trimming after adapter removal
  t5u="$sample_out/${base}_trim5_u${EXTRA_TRIM}x2.fastq"
  run_log "cutadapt -u $EXTRA_TRIM -u -$EXTRA_TRIM --minimum-length $MIN_LEN --cores $threads -o '$t5u' '$t5' > '$sample_out/${base}_cutadapt_u${EXTRA_TRIM}x2.log'"

  # alignment
  bam="$sample_out/${base}_sorted.bam"
  run_log "{ minimap2 -ax map-ont -t $threads '$mmi_index' '$t5u' | samtools view -Sb - | samtools sort -@ $threads -o '$bam'; }"
  run_log "samtools index '$bam'"

  # -------------------------
  # Variant calling (single callset)
  # -------------------------
  vcf_raw="$sample_out/${base}_all_variants_raw.vcf"
  vcf_all="$sample_out/${base}_all_variants.vcf"
  vcf_snps="$sample_out/${base}_snps.vcf"

  run_log "{ bcftools mpileup -d $PILEUP_MAX_DEPTH -Q$BASEQ_MIN -q$MAPQ_MIN -Ou -f '$ref' '$bam' \
    | bcftools call -cv --ploidy 1 -f GQ -Ov -o '$vcf_raw'; }"

  run_log "bcftools filter -i 'QUAL>'$QUAL_MIN'' -Ov -o '$vcf_all' '$vcf_raw'"

  # SNP-only output (derived from filtered callset)
  run_log "bcftools view -i 'TYPE=\"snp\"' -Ov -o '$vcf_snps' '$vcf_all'"

  # QC BAM
  run_log "fastqc '$bam' -o '$sample_out'"

  echo "Completed: $base" | tee -a "$log_file"
done

if ! $found_any; then
  echo "No FASTQ files found in '$input_folder'."
  exit 0
fi

echo "All outputs written to: $run_out"