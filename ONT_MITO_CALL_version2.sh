#!/opt/homebrew/bin/bash
# ONT_MITO_CALL_version2.sh
# tag: ONT_MITO_CALL_version2
#
# PURPOSE
# ---------------------------------------------------------------------------
# A specialized Nanopore mitochondrial DNA variant calling pipeline.
# It performs read quality filtering, rigorous trimming, alignment,
# and variant calling, followed by region-based annotation to produce
# four distinct VCF datasets for specific analytical needs.
#
# WORKFLOW OVERVIEW
# ---------------------------------------------------------------------------
# 1. READ FILTERING (QS Score)
#    - Scans FASTQ headers for 'qs:f:<val>' tags.
#    - Discards reads below QS_MIN (default: 10).
#    - Ensures high-confidence input data before processing begins.
#
# 2. PREPROCESSING
#    - FastQC (on QS-filtered reads).
#    - Cutadapt (Double-pass: 3' then 5' adapter trimming).
#    - Extra fixed trimming (EXTRA_TRIM) to remove potential edge artifacts.
#    - Length filtering (MIN_LEN).
#
# 3. ALIGNMENT
#    - Minimap2 (map-ont preset) against a linearized mitochondrial reference.
#    - Sorting and Indexing via Samtools.
#
# 4. VARIANT CALLING
#    - bcftools mpileup (High depth support: PILEUP_MAX_DEPTH).
#    - bcftools call (Ploidy 1).
#    - Quality filtering (QUAL > QUAL_MIN).
#
# 5. ANNOTATION & OUTPUT SPLITTING
#    Variants are annotated using a BED file to flag Homopolymer (HP)
#    regions and Blacklisted sites. The pipeline then generates:
#
#    A) *_annotated_all.vcf
#       - All called variants with region annotations.
#
#    B) *_snps.vcf
#       - SNPs only.
#       - EXCLUDES specific artifact positions (7898, 7899, 8595).
#
#    C) *_clean.vcf
#       - High-confidence variants.
#       - EXCLUDES Homopolymer regions and Blacklist sites.
#
#    D) *_homopolymers.vcf
#       - Contains ONLY variants found within Homopolymer regions.
#
# USAGE
# ---------------------------------------------------------------------------
#   ./ONT_MITO_CALL_version2.sh <input_fastq_folder>
#
# CONFIGURATION (Environment Variables)
# ---------------------------------------------------------------------------
#   export ref="linearized_mtdna.fasta"       # Reference genome
#   export regions_bed="linearized_regions.bed" # BED file for annotation
#   export QS_MIN=20                          # Min read quality score (qs:f)
#   export threads=8                          # CPU threads
#   export QUAL_MIN=20                        # Variant QUAL score cutoff
#   export MIN_LEN=90                         # Min read length
#   export EXTRA_TRIM=22                      # Fixed trim from ends
#   export PILEUP_MAX_DEPTH=100000            # Max depth for pileup
#   export BASEQ_MIN=20                       # Min base quality
#   export MAPQ_MIN=20                        # Min mapping quality
# ---------------------------------------------------------------------------

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <input_fastq_folder>"
  exit 1
fi

input_folder="$1"
[[ -d "$input_folder" ]] || { echo "Error: Input folder '$input_folder' not found."; exit 1; }

# --- Configuration ---
threads="${threads:-8}"
pipeline_name="ONT_MITO_CALL_version2"

# Keep reads with qs >= QS_MIN (filter out reads below)
QS_MIN="${QS_MIN:-10}"

ref="${ref:-linearized_mtdna.fasta}"
mmi_index="${mmi_index:-${ref}.mmi}"
regions_bed="${regions_bed:-linearized_regions.bed}"

QUAL_MIN="${QUAL_MIN:-20}"
MIN_LEN="${MIN_LEN:-90}"
EXTRA_TRIM="${EXTRA_TRIM:-22}"
PILEUP_MAX_DEPTH="${PILEUP_MAX_DEPTH:-100000}"
BASEQ_MIN="${BASEQ_MIN:-20}"
MAPQ_MIN="${MAPQ_MIN:-20}"

# --- Checks ---
if [[ ! -f "$regions_bed" ]]; then
  echo "Error: '$regions_bed' not found."
  echo "Please run the python script first to generate the regions BED file."
  exit 1
fi

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
  echo "Error: adapter file not found"
  echo "Looked for: Updated_Adapter_Primer_List_Cutadapt_cleaned.txt"
  exit 1
fi

run_ts="$(date +%Y%m%d_%H%M%S)"
run_out="${pipeline_name}_${run_ts}_output"
mkdir -p "$run_out"

# --- Get Reference Chromosome Name ---
ref_chrom_name=$(head -n1 "$ref" | cut -d ' ' -f1 | tr -d '>')
echo "Detected Reference Chromosome Name: $ref_chrom_name"

run_log() {
  # Run external commands with logging, no nested bash -lc needed
  {
    printf 'Running: %s\n' "$*"
    "$@"
    printf '\n'
  } 2>&1 | tee -a "$log_file"
}

# FASTQ filter: keep reads with qs:f:<val> >= QS_MIN
# Works with .fastq or .fastq.gz; assumes 4-line FASTQ.
filter_qs_fastq() {
  local in_fq="$1"
  local out_fq="$2"
  local qs_min="$3"

  local reader=(cat)
  [[ "$in_fq" == *.gz ]] && reader=(gzip -cd)

  "${reader[@]}" "$in_fq" | awk -v qs_min="$qs_min" '
    function get_qs(h,   n,i,a,v) {
      n = split(h, a, "\t");
      for (i=1; i<=n; i++) {
        if (index(a[i], "qs:f:") == 1) {
          v = substr(a[i], 6);
          return v + 0;
        }
      }
      return -1;
    }
    NR%4==1 {
      h=$0;
      qs=get_qs(h);
      keep = (qs >= qs_min);
      if (keep) print h;
      next
    }
    NR%4==2 { if (keep) print $0; next }
    NR%4==3 { if (keep) print $0; next }
    NR%4==0 { if (keep) print $0; next }
  ' > "$out_fq"
}

# Reference indexing
[[ -f "$mmi_index" ]] || minimap2 -d "$mmi_index" "$ref"
[[ -f "${ref}.fai" ]] || samtools faidx "$ref"

shopt -s nullglob
found_any=false

for fq in "$input_folder"/*.fastq "$input_folder"/*.fastq.gz; do
  [[ -f "$fq" ]] || continue
  found_any=true

  base="$(basename "$fq")"
  base="${base%.fastq}"
  base="${base%.fastq.gz}"

  sample_out="$run_out/${base}"
  mkdir -p "$sample_out"
  log_file="$sample_out/${base}.log"

  echo "Processing $base" | tee -a "$log_file"
  echo "QS filter: KEEP reads with qs >= $QS_MIN (FILTER OUT qs < $QS_MIN)" | tee -a "$log_file"

  # --- QS filtering (run directly, log separately) ---
  qs_fq="$sample_out/${base}_qsGE${QS_MIN}.fastq"
  {
    echo "Running: filter_qs_fastq '$fq' '$qs_fq' '$QS_MIN'"
    filter_qs_fastq "$fq" "$qs_fq" "$QS_MIN"
    echo
  } 2>&1 | tee -a "$log_file"

  # Count kept reads
  run_log awk 'END{print "Kept reads (approx):", NR/4}' "$qs_fq"

  # --- Preprocessing ---
  run_log fastqc "$qs_fq" -o "$sample_out"

  t3="$sample_out/${base}_qsGE${QS_MIN}_trim3.fastq"
  t5="$sample_out/${base}_qsGE${QS_MIN}_trim5.fastq"
  run_log cutadapt -a "file:$ADAPTER_FILE" --error-rate 0.10 --overlap 5 --minimum-length "$MIN_LEN" --cores "$threads" -o "$t3" "$qs_fq"
  run_log bash -lc "echo '(cutadapt -a log) written by cutadapt above'"

  run_log cutadapt -g "file:$ADAPTER_FILE" --error-rate 0.10 --overlap 5 --minimum-length "$MIN_LEN" --cores "$threads" -o "$t5" "$t3"

  t5u="$sample_out/${base}_qsGE${QS_MIN}_trim5_u${EXTRA_TRIM}x2.fastq"
  run_log cutadapt -u "$EXTRA_TRIM" -u "-$EXTRA_TRIM" --minimum-length "$MIN_LEN" --cores "$threads" -o "$t5u" "$t5"

  # --- Alignment ---
  bam="$sample_out/${base}_qsGE${QS_MIN}_sorted.bam"
  run_log bash -lc "{ minimap2 -ax map-ont -t $threads '$mmi_index' '$t5u' | samtools view -Sb - | samtools sort -@ $threads -o '$bam'; }"
  run_log samtools index "$bam"

  # --- Variant Calling ---
  vcf_raw="$sample_out/${base}_qsGE${QS_MIN}_raw.vcf"
  vcf_qual="$sample_out/${base}_qsGE${QS_MIN}_qual_filtered.vcf"

  run_log bash -lc "{ bcftools mpileup -d $PILEUP_MAX_DEPTH -Q$BASEQ_MIN -q$MAPQ_MIN -Ou -f '$ref' '$bam' | bcftools call -cv --ploidy 1 -f GQ -Ov -o '$vcf_raw'; }"
  run_log bcftools filter -i "QUAL>$QUAL_MIN" -Ov -o "$vcf_qual" "$vcf_raw"

  # --- Annotation & Subsetting ---
  sample_bed="$sample_out/${base}_targets.bed"
  sed "s/PLACEHOLDER_CHROM/$ref_chrom_name/g" "$regions_bed" > "$sample_bed"

  echo '##INFO=<ID=RegionType,Number=1,Type=String,Description="Type of region (HP_Region or Blacklist_Site)">' > "$sample_out/hdr.txt"
  vcf_annotated="$sample_out/${base}_qsGE${QS_MIN}_annotated_all.vcf"
  run_log bcftools annotate -a "$sample_bed" -c CHROM,FROM,TO,RegionType -h "$sample_out/hdr.txt" -Ov -o "$vcf_annotated" "$vcf_qual"

  vcf_snps="$sample_out/${base}_qsGE${QS_MIN}_snps.vcf"
  run_log bash -lc "bcftools view -i 'TYPE=\"snp\"' '$vcf_annotated' | bcftools filter -e 'POS=7898 || POS=7899 || POS=8595' -Ov -o '$vcf_snps'"

  vcf_clean="$sample_out/${base}_qsGE${QS_MIN}_clean.vcf"
  run_log bcftools filter -e 'RegionType="HP_Region" || RegionType="Blacklist_Site"' -Ov -o "$vcf_clean" "$vcf_annotated"

  vcf_hp="$sample_out/${base}_qsGE${QS_MIN}_homopolymers.vcf"
  run_log bcftools filter -i 'RegionType="HP_Region"' -Ov -o "$vcf_hp" "$vcf_annotated"

  rm -f "$vcf_raw" "$vcf_qual" "$sample_bed" "$sample_out/hdr.txt"

  # QC (keeping your original behavior)
  run_log fastqc "$bam" -o "$sample_out"

  echo "Completed: $base" | tee -a "$log_file"
  echo "  -> QS-filtered FASTQ: $qs_fq" | tee -a "$log_file"
  echo "  -> Annotated (All):   $vcf_annotated" | tee -a "$log_file"
  echo "  -> SNPs (Special):    $vcf_snps" | tee -a "$log_file"
  echo "  -> Clean (No HP):     $vcf_clean" | tee -a "$log_file"
  echo "  -> Homopolymers:      $vcf_hp" | tee -a "$log_file"
done

if ! $found_any; then
  echo "No FASTQ files found in '$input_folder' (expected *.fastq or *.fastq.gz)."
  exit 0
fi

echo "All outputs written to: $run_out"