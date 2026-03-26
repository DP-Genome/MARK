#!/opt/homebrew/bin/bash
# ILM_MITO_v32.sh
# tag: ILM_MITO_v32
#
# PURPOSE
# ---------------------------------------------------------------------------
# Illumina Merged-Read Diagnostic Split-Track Pipeline for mitochondrial data.
# 
# CHANGES in v32:
# - Retains the original midpoint-assignment masking logic from v31.
# - Adds a minimal EDGE_PROTECT window at each assigned amplicon boundary to preserve true overlap-edge variants.
# - Keeps the original masking model everywhere else to avoid broad callset disruption.
# - Removed FASTA reference loading from Python memory to speed up processing.
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
pipeline_name="ILM_MITO_v32"

# LENSAFE setting (Default 3000 for Illumina)
MAX_LEN="${MAX_LEN:-3000}"
DISCARD_WARN_PCT="${DISCARD_WARN_PCT:-5}"

# fastp settings
READQ="${READQ:-20}"
UNQUAL_PCT="${UNQUAL_PCT:-40}"
MIN_LEN="${MIN_LEN:-90}"
N_BASE_LIMIT="${N_BASE_LIMIT:-5}"

# cutadapt settings
CUTADAPT_ERR="${CUTADAPT_ERR:-0.10}"
CUTADAPT_OVL="${CUTADAPT_OVL:-5}"
EXTRA_TRIM="${EXTRA_TRIM:-0}"

ref="${ref:-linearized_mtdna.fasta}"
regions_bed="${regions_bed:-linearized_regions.bed}"

# HIERARCHICAL FILTERS
SAFETY_QUAL="${SAFETY_QUAL:-20}"
STRICT_QUAL="${STRICT_QUAL:-60}"
MIN_DEPTH="${MIN_DEPTH:-10}"       
PILEUP_MAX_DEPTH="${PILEUP_MAX_DEPTH:-100000}"
BASEQ_MIN="${BASEQ_MIN:-20}"       
MAPQ_MIN="${MAPQ_MIN:-20}"
EDGE_PROTECT="${EDGE_PROTECT:-1}"

# --- Checks ---
if [[ ! -f "$regions_bed" ]]; then
  echo "Error: '$regions_bed' not found."
  exit 1
fi
if [[ ! -f "$ref" ]]; then
  echo "Error: Reference FASTA '$ref' not found."
  exit 1
fi

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
  exit 1
fi

# Ensure reference indices exist
need_bwa_index=false
for ext in .0123 .amb .ann .bwt.2bit.64 .pac; do
  if [[ ! -f "${ref}${ext}" ]]; then
    need_bwa_index=true
    break
  fi
done
if $need_bwa_index; then
  echo "Building bwa-mem2 index for: $ref"
  bwa-mem2 index "$ref"
fi
[[ -f "${ref}.fai" ]] || samtools faidx "$ref"

command -v fastp >/dev/null 2>&1 || { echo "Error: fastp not found."; exit 1; }

# --- Output Folder Setup ---
run_ts="$(date +%Y%m%d_%H%M%S)"
input_name="$(basename "$input_folder")"
run_out="${pipeline_name}_${input_name}_${run_ts}_output"
mkdir -p "$run_out"

ref_chrom_name=$(head -n1 "$ref" | cut -d ' ' -f1 | tr -d '>')

run_log() {
  {
    printf 'Running: %s\n' "$*"
    "$@"
    printf '\n'
  } 2>&1 | tee -a "$log_file"
}

count_fastq_reads() {
  local fq="$1"
  if [[ "$fq" == *.gz ]]; then
    gzip -cd "$fq" | awk 'END{printf "%.0f\n", NR/4}'
  else
    awk 'END{printf "%.0f\n", NR/4}' "$fq"
  fi
}

shopt -s nullglob
found_any=false

# Global Run summary
run_summary_file="$run_out/run_summary.txt"
{
  echo "## =============================================================================="
  echo "## PIPELINE RUN SUMMARY: DIAGNOSTIC RETENTION TRACKING"
  echo "## =============================================================================="
  echo "## Date Run       : $(date)"
  echo "## Pipeline Name  : $pipeline_name"
  echo "## Input Folder   : $(realpath "$input_folder")"
  echo "## Threads        : $threads"
  echo "## Reference      : $(realpath "$ref")"
  echo "## Regions BED    : $(realpath "$regions_bed")"
  echo "## Adapter File   : $(realpath "$ADAPTER_FILE")"
  echo "##"
  echo "## --- PRE-PROCESSING & TRIMMING ---"
  echo "## MAX_LEN        : $MAX_LEN"
  echo "## FASTP MERGE    : READQ=$READQ | UNQUAL_PCT=$UNQUAL_PCT | MIN_LEN=$MIN_LEN"
  echo "## MIXED MATCH    : ERR=$CUTADAPT_ERR | OVL=$CUTADAPT_OVL | MIN_LEN=$MIN_LEN "
  echo "## EXTRA_TRIM     : $EXTRA_TRIM"
  echo "##"
  echo "## --- ALIGNMENT & CALLING ---"
  echo "## MAPQ_MIN       : $MAPQ_MIN"
  echo "## BASEQ_MIN      : $BASEQ_MIN"
  echo "## PILEUP_DEPTH   : $PILEUP_MAX_DEPTH"
  echo "##"
  echo "## --- VARIANT FILTERING ---"
  echo "## MIN_DEPTH      : $MIN_DEPTH"
  echo "## SAFETY_QUAL    : $SAFETY_QUAL"
  echo "## STRICT_QUAL    : $STRICT_QUAL"
  echo "##"
  echo "## =============================================================================="
  printf "sample\tstage\treads_in\treads_out\tdropped_reads\tdropped_pct\n"
} > "$run_summary_file"

for r1 in \
  "$input_folder"/*_R1*.fastq "$input_folder"/*_R1*.fastq.gz \
  "$input_folder"/*.R1*.fastq "$input_folder"/*.R1*.fastq.gz
do
  [[ -f "$r1" ]] || continue

  r2="$r1"
  r2="${r2/_R1/_R2}"
  r2="${r2/.R1/.R2}"

  if [[ ! -f "$r2" ]]; then
    echo "Warning: Missing R2 pair for R1: $r1"
    continue
  fi

  found_any=true

  base="$(basename "$r1")"
  base="${base%.fastq}"
  base="${base%.fastq.gz}"
  base="${base%_R1*}"
  base="${base%.R1*}"

  sample_out="$run_out/${base}"
  mkdir -p "$sample_out"
  log_file="$sample_out/${base}.log"

  echo "Processing $base" | tee -a "$log_file"
  input_reads=$(count_fastq_reads "$r1") # Multiply by 2 for total reads later if wanted, tracking pairs here

  # --- STEP 0: LENSAFE ---
  lensafe_r1="$sample_out/${base}_lenLE${MAX_LEN}_R1.fastq"
  lensafe_r2="$sample_out/${base}_lenLE${MAX_LEN}_R2.fastq"
  fp_lensafe_html="$sample_out/${base}_fastp_lensafe.html"
  fp_lensafe_json="$sample_out/${base}_fastp_lensafe.json"
  
  run_log fastp \
    -i "$r1" -I "$r2" \
    -o "$lensafe_r1" -O "$lensafe_r2" \
    --length_limit "$MAX_LEN" \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --thread "$threads" \
    --html "$fp_lensafe_html" \
    --json "$fp_lensafe_json" 2> /dev/null
    
  len_reads=$(count_fastq_reads "$lensafe_r1")
  len_drop=$((input_reads - len_reads))
  len_pct=$(awk -v d="$len_drop" -v i="$input_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t0_Lensafe(Pairs)\t%s\t%s\t%s\t%s%%\n" "$base" "$input_reads" "$len_reads" "$len_drop" "$len_pct" >> "$run_summary_file"

  # --- STEP 1: Merge & Filter ---
  merged_fq="$sample_out/${base}_merged.fastq"
  fp_html="$sample_out/${base}_fastp_merge.html"
  fp_json="$sample_out/${base}_fastp_merge.json"

  run_log fastp \
    -i "$lensafe_r1" -I "$lensafe_r2" \
    --merge --include_unmerged \
    --merged_out "$merged_fq" \
    --qualified_quality_phred "$READQ" \
    --unqualified_percent_limit "$UNQUAL_PCT" \
    --length_required "$MIN_LEN" \
    --n_base_limit "$N_BASE_LIMIT" \
    --thread "$threads" \
    --html "$fp_html" \
    --json "$fp_json" 2> /dev/null

  merged_reads=$(count_fastq_reads "$merged_fq")
  merged_drop=$((len_reads * 2 - merged_reads))
  merged_pct=$(awk -v d="$merged_drop" -v i="$((len_reads * 2))" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t1_MergeFilter(Reads)\t%s\t%s\t%s\t%s%%\n" "$base" "$((len_reads * 2))" "$merged_reads" "$merged_drop" "$merged_pct" >> "$run_summary_file"

  # --- STEP 2: Preprocessing ---
  run_log fastqc "$merged_fq" -o "$sample_out" > /dev/null 2>&1
  
  t3="$sample_out/${base}_trim3.fastq"
  run_log cutadapt -a "file:$ADAPTER_FILE" --error-rate "$CUTADAPT_ERR" --overlap "$CUTADAPT_OVL" --minimum-length "$MIN_LEN" --cores "$threads" -o "$t3" "$merged_fq" >/dev/null
  t3_reads=$(count_fastq_reads "$t3")
  t3_drop=$((merged_reads - t3_reads))
  t3_pct=$(awk -v d="$t3_drop" -v i="$merged_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t2_Trim3\t%s\t%s\t%s\t%s%%\n" "$base" "$merged_reads" "$t3_reads" "$t3_drop" "$t3_pct" >> "$run_summary_file"

  t5="$sample_out/${base}_trim5.fastq"
  run_log cutadapt -g "file:$ADAPTER_FILE" --error-rate "$CUTADAPT_ERR" --overlap "$CUTADAPT_OVL" --minimum-length "$MIN_LEN" --cores "$threads" -o "$t5" "$t3" >/dev/null
  t5_reads=$(count_fastq_reads "$t5")
  t5_drop=$((t3_reads - t5_reads))
  t5_pct=$(awk -v d="$t5_drop" -v i="$t3_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t3_Trim5\t%s\t%s\t%s\t%s%%\n" "$base" "$t3_reads" "$t5_reads" "$t5_drop" "$t5_pct" >> "$run_summary_file"
  
  final_fq="$sample_out/${base}_trim5_u${EXTRA_TRIM}x2.fastq"
  if [[ "$EXTRA_TRIM" -gt 0 ]]; then
    run_log cutadapt -u "$EXTRA_TRIM" -u "-$EXTRA_TRIM" --minimum-length "$MIN_LEN" --cores "$threads" -o "$final_fq" "$t5" >/dev/null
  else
    ln -sf "$(basename "$t5")" "$final_fq"
  fi
  final_reads=$(count_fastq_reads "$final_fq")
  final_drop=$((t5_reads - final_reads))
  final_pct=$(awk -v d="$final_drop" -v i="$t5_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t4_ExtraTrim\t%s\t%s\t%s\t%s%%\n" "$base" "$t5_reads" "$final_reads" "$final_drop" "$final_pct" >> "$run_summary_file"


  # --- STEP 3: Initial Alignment ---
  bam_initial="$sample_out/${base}_initial_sorted.bam"
  run_log bash -lc "{ bwa-mem2 mem -t $threads '$ref' '$final_fq' | samtools view -Sb - | samtools sort -@ $threads -o '$bam_initial'; }"
  run_log samtools index "$bam_initial"

  # =========================================================================
  # TRACK A: BASELINE (NO MASKING)
  # =========================================================================
  echo "[Pipeline] Generating BASELINE track..." | tee -a "$log_file"
  bam_baseline="$sample_out/${base}_baseline_sorted.bam"
  cp "$bam_initial" "$bam_baseline"
  cp "${bam_initial}.bai" "${bam_baseline}.bai"

  # Baseline Variant Calling
  vcf_raw_base="$sample_out/${base}_baseline_raw.vcf"
  vcf_qual_base="$sample_out/${base}_baseline_qual_filtered.vcf"
  run_log bash -lc "{ bcftools mpileup -a FORMAT/AD,FORMAT/DP -d $PILEUP_MAX_DEPTH -Q$BASEQ_MIN -q$MAPQ_MIN -Ou -f '$ref' '$bam_baseline' | bcftools call -mv --ploidy 1 -Ov -o '$vcf_raw_base'; }"
  run_log bcftools filter -i "QUAL>$SAFETY_QUAL && INFO/DP>=$MIN_DEPTH" -Ov -o "$vcf_qual_base" "$vcf_raw_base"

  # Baseline Annotation
  if [[ -s "$vcf_qual_base" ]]; then
    sample_bed="$sample_out/${base}_targets.bed"
    sed "s/PLACEHOLDER_CHROM/$ref_chrom_name/g" "$regions_bed" > "$sample_bed"
    echo '##INFO=<ID=RegionType,Number=1,Type=String,Description="Type of region (HP_Region or Blacklist_Site)">' > "$sample_out/hdr.txt"

    vcf_annotated_base="$sample_out/${base}_baseline_annotated_all.vcf"
    run_log bcftools annotate -a "$sample_bed" -c CHROM,FROM,TO,RegionType -h "$sample_out/hdr.txt" -Ov -o "$vcf_annotated_base" "$vcf_qual_base"

    run_log bash -lc "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_base' | bcftools filter -e 'POS=7898 || POS=7899 || POS=8595' -Ov -o '$sample_out/${base}_baseline_snps.vcf'"
    run_log bash -lc "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_base' | bcftools filter -e 'RegionType=\"HP_Region\" || RegionType=\"Blacklist_Site\"' -Ov -o '$sample_out/${base}_baseline_clean.vcf'"
    run_log bcftools filter -i 'RegionType="HP_Region"' -Ov -o "$sample_out/${base}_baseline_homopolymers.vcf" "$vcf_annotated_base"
  fi


  # =========================================================================
  # TRACK B: EXPERIMENTAL (PYTHON MASKING)
  # =========================================================================
  echo "[Pipeline] Generating MASKED track..." | tee -a "$log_file"
  bam_masked="$sample_out/${base}_masked_sorted.bam"

  samtools view -h "$bam_initial" | python3 -c '
import os, sys, re

# v32 overlap-rescue boundary protection using original midpoint-assignment amplicons (Linearized Coordinates)
amplicons = [
    ("Amp1", 7729, 7842),
    ("Amp2", 7832, 7941),
    ("Amp3", 7939, 8124),
    ("Amp4", 8103, 8202),
    ("Amp5", 8190, 8315),
    ("Amp6", 8271, 8437),
    ("Amp7", 8421, 8542),
    ("Amp8", 8531, 8649),
    ("Amp9", 8627, 8721),
    ("Amp10", 8714, 8877)
]

cigar_re = re.compile(r"(\d+)([MIDNSHPX=])")
edge_protect = int(os.environ.get("EDGE_PROTECT", "1"))

for line in sys.stdin:
    if line.startswith("@"):
        sys.stdout.write(line)
        continue

    line = line.rstrip("\n")
    parts = line.split("\t")
    if len(parts) < 11: 
        sys.stdout.write(line + "\n")
        continue

    flag = int(parts[1])
    # Skip unmapped and secondary alignments to prevent leaks
    if flag & 4 or flag & 256:
        sys.stdout.write(line + "\n")
        continue

    pos = int(parts[3]) - 1 # 0-based
    cigar = parts[5]
    seq = parts[9]
    qual = parts[10]
    
    if cigar == "*" or qual == "*" or seq == "*":
        sys.stdout.write(line + "\n")
        continue

    ref_span = sum(int(n) for n, op in cigar_re.findall(cigar) if op in "MDN=X")
    ref_end = pos + ref_span
    
    # --- V32: ORIGINAL MIDPOINT ASSIGNMENT + MINIMAL EDGE PROTECTION ---
    read_midpoint = (pos + ref_end) / 2
    best_amp = None
    min_dist = float("inf")
    
    for name, a_start, a_end in amplicons:
        amp_midpoint = (a_start + a_end) / 2
        dist = abs(read_midpoint - amp_midpoint)
        if dist < min_dist:
            min_dist = dist
            best_amp = (a_start, a_end)
            
    target_start, target_end = best_amp
    keep_start = target_start - edge_protect
    keep_end = target_end + edge_protect

    if ref_end < keep_start or pos > keep_end:
        continue 

    ref_p = pos
    read_p = 0
    qual_list = list(qual)
    
    for n_str, op in cigar_re.findall(cigar):
        n = int(n_str)
        if op in "M=X":
            for _ in range(n):
                if not (keep_start <= ref_p < keep_end):
                    if read_p < len(qual_list):
                        qual_list[read_p] = "!"
                ref_p += 1; read_p += 1
        elif op in "IS":
            read_p += n
        elif op in "DN":
            ref_p += n
            
    parts[10] = "".join(qual_list)
    try:
        sys.stdout.write("\t".join(parts) + "\n")
    except BrokenPipeError:
        sys.exit(0)
' | samtools view -Sb - | samtools sort -@ "$threads" -o "$bam_masked"

  run_log samtools index "$bam_masked"

  # Masked Variant Calling
  vcf_raw_masked="$sample_out/${base}_masked_raw.vcf"
  vcf_qual_masked="$sample_out/${base}_masked_qual_filtered.vcf"
  run_log bash -lc "{ bcftools mpileup -a FORMAT/AD,FORMAT/DP -d $PILEUP_MAX_DEPTH -Q$BASEQ_MIN -q$MAPQ_MIN -Ou -f '$ref' '$bam_masked' | bcftools call -mv --ploidy 1 -Ov -o '$vcf_raw_masked'; }"
  run_log bcftools filter -i "QUAL>$SAFETY_QUAL && INFO/DP>=$MIN_DEPTH" -Ov -o "$vcf_qual_masked" "$vcf_raw_masked"

  # Masked Annotation
  if [[ -s "$vcf_qual_masked" ]]; then
    vcf_annotated_masked="$sample_out/${base}_masked_annotated_all.vcf"
    run_log bcftools annotate -a "$sample_bed" -c CHROM,FROM,TO,RegionType -h "$sample_out/hdr.txt" -Ov -o "$vcf_annotated_masked" "$vcf_qual_masked"

    run_log bash -lc "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_masked' | bcftools filter -e 'POS=7898 || POS=7899 || POS=8595' -Ov -o '$sample_out/${base}_masked_snps.vcf'"
    run_log bash -lc "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_masked' | bcftools filter -e 'RegionType=\"HP_Region\" || RegionType=\"Blacklist_Site\"' -Ov -o '$sample_out/${base}_masked_clean.vcf'"
    run_log bcftools filter -i 'RegionType="HP_Region"' -Ov -o "$sample_out/${base}_masked_homopolymers.vcf" "$vcf_annotated_masked"
    
    rm -f "$sample_bed" "$sample_out/hdr.txt"
  fi

  # =========================================================================
  # DIAGNOSTIC DEPTH COMPARISON
  # =========================================================================
  echo "[Diagnostic] Calculating average depth across genome..." | tee -a "$log_file"
  depth_baseline=$(samtools depth -a -q "$BASEQ_MIN" "$bam_baseline" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
  depth_masked=$(samtools depth -a -q "$BASEQ_MIN" "$bam_masked" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
  
  echo "---------------------------------------------------" | tee -a "$log_file"
  echo "DIAGNOSTIC SUMMARY: $base" | tee -a "$log_file"
  echo "Baseline Average Depth : $depth_baseline" | tee -a "$log_file"
  echo "Masked Average Depth   : $depth_masked" | tee -a "$log_file"
  echo "---------------------------------------------------" | tee -a "$log_file"

  run_log fastqc "$bam_masked" -o "$sample_out" > /dev/null 2>&1
  echo "Completed Diagnostic Split: $base" | tee -a "$log_file"
done

if ! $found_any; then
  echo "No paired FASTQs found in '$input_folder'."
  exit 0
fi

echo "All outputs written to: $run_out"