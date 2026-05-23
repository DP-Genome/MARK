#!/usr/bin/env bash
# MARK.sh
# tag: MARK
#
# PURPOSE
# ---------------------------------------------------------------------------
# Diagnostic Split-Track Pipeline for Nanopore mitochondrial data.
# 
# CHANGES in v36:
# - Fixed 1bp gap on the left border of amplicons (0-based conversion constraint).
# - Retains v34/v35 CIGAR soft-clipping and maximum overlap logic.
# ---------------------------------------------------------------------------

set -euo pipefail

VERSION="1.1"

if [[ $# -lt 1 || "$1" == "-h" || "$1" == "--help" ]]; then
  echo "MARK Pipeline v$VERSION"
  echo "Usage: $0 <input_fastq_folder>"
  exit 1
fi

input_folder="$1"
[[ -d "$input_folder" ]] || { echo "Error: Input folder '$input_folder' not found."; exit 1; }

# --- Configuration ---
threads="${threads:-8}"
pipeline_name="MARK"

MAX_LEN="${MAX_LEN:-1500}"
DISCARD_WARN_PCT="${DISCARD_WARN_PCT:-5}"
QS_MIN="${QS_MIN:-10}"
MIN_LEN="${MIN_LEN:-90}"
MIN_LEN_POST="${MIN_LEN_POST:-90}"
MAX_LEN_POST="${MAX_LEN_POST:-300}"
EXTRA_TRIM="${EXTRA_TRIM:-0}"
CUTADAPT_ERR="${CUTADAPT_ERR:-0.10}" 
CUTADAPT_OVL="${CUTADAPT_OVL:-5}"    

ref="${ref:-linearized_mtdna.fasta}"
mmi_index="${mmi_index:-${ref}.mmi}"
regions_bed="${regions_bed:-linearized_regions.bed}"

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

command -v fastplong >/dev/null 2>&1 || { echo "Error: fastplong not found."; exit 1; }

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

detect_fastplong_maxlen_flag() {
  local help
  help="$(fastplong --help 2>&1 || true)"
  local candidates=("--max_len" "--max_length" "--max_read_length" "--length_limit" "--length_max" "--max_length_required")
  for f in "${candidates[@]}"; do
    if echo "$help" | grep -qE "^\s*${f}\b|[[:space:]]${f}\b"; then echo "$f"; return 0; fi
  done
  return 1
}

FASTPLONG_MAXLEN_FLAG="$(detect_fastplong_maxlen_flag || true)"

[[ -f "$mmi_index" ]] || minimap2 -d "$mmi_index" "$ref"
[[ -f "${ref}.fai" ]] || samtools faidx "$ref"

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
  echo "## QS_MIN         : $QS_MIN"
  echo "## MIXED MATCH    : ERR=$CUTADAPT_ERR | OVL=$CUTADAPT_OVL | MIN_LEN=$MIN_LEN "
  echo "## EXTRA_TRIM     : $EXTRA_TRIM"
  echo "## MAX_LEN_POST   : $MAX_LEN_POST"
  echo "## MIN_LEN_POST   : $MIN_LEN_POST"
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

for fq in "$input_folder"/*.fastq "$input_folder"/*.fastq.gz; do
  [[ -f "$fq" ]] || continue
  found_any=true

  base="$(basename "$fq" | sed 's/\.fastq.*//')"
  sample_out="$run_out/${base}"
  mkdir -p "$sample_out"
  log_file="$sample_out/${base}.log"

  echo "Processing $base" | tee -a "$log_file"
  input_reads=$(count_fastq_reads "$fq")

  # --- STEP 0: LENSAFE ---
  maxlen_fq="$sample_out/${base}_lenLE${MAX_LEN}.fastq"
  maxlen_failed="$sample_out/${base}_lenGT${MAX_LEN}_FAILED.fastq"
  fp_html="$sample_out/${base}_fastplong.html"
  fp_json="$sample_out/${base}_fastplong.json"
  
  run_log fastplong -i "$fq" -o "$maxlen_fq" --failed_out "$maxlen_failed" "$FASTPLONG_MAXLEN_FLAG" "$MAX_LEN" -h "$fp_html" -j "$fp_json"
  len_reads=$(count_fastq_reads "$maxlen_fq")
  len_drop=$((input_reads - len_reads))
  len_pct=$(awk -v d="$len_drop" -v i="$input_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t0_Lensafe\t%s\t%s\t%s\t%s%%\n" "$base" "$input_reads" "$len_reads" "$len_drop" "$len_pct" >> "$run_summary_file"

  # --- STEP 1: QS Filtering ---
  qs_fq="$sample_out/${base}_qsGE${QS_MIN}.fastq"
  if head -n 1 "$maxlen_fq" | grep -q "qs:f:"; then
      awk -v qs_min="$QS_MIN" '
        function get_qs(h, n,i,a,v) {
          n = split(h, a, "\t");
          for (i=1; i<=n; i++) { if (index(a[i], "qs:f:") == 1) return substr(a[i], 6) + 0; }
          return -1;
        }
        NR%4==1 { h=$0; keep = (get_qs(h) >= qs_min); if(keep) print h; next }
        { if(keep) print $0 }
      ' "$maxlen_fq" > "$qs_fq"
  else
      cp "$maxlen_fq" "$qs_fq"
  fi
  qs_reads=$(count_fastq_reads "$qs_fq")
  qs_drop=$((len_reads - qs_reads))
  qs_pct=$(awk -v d="$qs_drop" -v i="$len_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t1_QS_GE%s\t%s\t%s\t%s\t%s%%\n" "$base" "$QS_MIN" "$len_reads" "$qs_reads" "$qs_drop" "$qs_pct" >> "$run_summary_file"

  # --- STEP 2: Preprocessing ---
  run_log fastqc "$qs_fq" -o "$sample_out"
  
  t3="$sample_out/${base}_qsGE${QS_MIN}_trim3.fastq"
  run_log cutadapt -a "file:$ADAPTER_FILE" --error-rate "$CUTADAPT_ERR" --overlap "$CUTADAPT_OVL" --minimum-length "$MIN_LEN" --cores "$threads" -o "$t3" "$qs_fq" >/dev/null
  t3_reads=$(count_fastq_reads "$t3")
  t3_drop=$((qs_reads - t3_reads))
  t3_pct=$(awk -v d="$t3_drop" -v i="$qs_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t2_Trim3\t%s\t%s\t%s\t%s%%\n" "$base" "$qs_reads" "$t3_reads" "$t3_drop" "$t3_pct" >> "$run_summary_file"

  t5="$sample_out/${base}_qsGE${QS_MIN}_trim5.fastq"
  run_log cutadapt -g "file:$ADAPTER_FILE" --error-rate "$CUTADAPT_ERR" --overlap "$CUTADAPT_OVL" --minimum-length "$MIN_LEN" --cores "$threads" -o "$t5" "$t3" >/dev/null
  t5_reads=$(count_fastq_reads "$t5")
  t5_drop=$((t3_reads - t5_reads))
  t5_pct=$(awk -v d="$t5_drop" -v i="$t3_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t3_Trim5\t%s\t%s\t%s\t%s%%\n" "$base" "$t3_reads" "$t5_reads" "$t5_drop" "$t5_pct" >> "$run_summary_file"
  
  final_fq="$sample_out/${base}_qsGE${QS_MIN}_trim5_u${EXTRA_TRIM}x2.fastq"
  if [[ "$EXTRA_TRIM" -gt 0 ]]; then
    run_log cutadapt -u "$EXTRA_TRIM" -u "-$EXTRA_TRIM" --minimum-length "$MIN_LEN" --cores "$threads" -o "$final_fq" "$t5" >/dev/null
  else
    ln -sf "$(basename "$t5")" "$final_fq"
  fi
  final_reads=$(count_fastq_reads "$final_fq")
  final_drop=$((t5_reads - final_reads))
  final_pct=$(awk -v d="$final_drop" -v i="$t5_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t4_ExtraTrim\t%s\t%s\t%s\t%s%%\n" "$base" "$t5_reads" "$final_reads" "$final_drop" "$final_pct" >> "$run_summary_file"

  # --- STEP 2.5: Post-Trim Length Filtering ---
  final_len_fq="$sample_out/${base}_qsGE${QS_MIN}_trim5_u${EXTRA_TRIM}x2_lenFiltered.fastq"
  run_log fastp -i "$final_fq" -o "$final_len_fq" \
    --length_required "$MIN_LEN_POST" \
    --length_limit "$MAX_LEN_POST" \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --thread "$threads" \
    -h "$sample_out/${base}_fastp_post_len.html" \
    -j "$sample_out/${base}_fastp_post_len.json" 2> /dev/null

  lenpost_reads=$(count_fastq_reads "$final_len_fq")
  lenpost_drop=$((final_reads - lenpost_reads))
  lenpost_pct=$(awk -v d="$lenpost_drop" -v i="$final_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t5_PostTrim_Len\t%s\t%s\t%s\t%s%%\n" "$base" "$final_reads" "$lenpost_reads" "$lenpost_drop" "$lenpost_pct" >> "$run_summary_file"
  
  final_fq="$final_len_fq"


  # --- STEP 3: Initial Alignment ---
  bam_initial="$sample_out/${base}_initial_sorted.bam"
  run_log bash -lc "{ minimap2 -ax map-ont -t $threads '$mmi_index' '$final_fq' | samtools view -Sb - | samtools sort -@ $threads -o '$bam_initial'; }"
  run_log samtools index "$bam_initial"

  # =========================================================================
  # TRACK A: BASELINE (NO TRIMMING)
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
  # TRACK B: EXPERIMENTAL (PYTHON TRIMMING/SOFT-CLIPPING)
  # =========================================================================
  echo "[Pipeline] Generating TRIMMED track..." | tee -a "$log_file"
  bam_trimmed="$sample_out/${base}_trimmed_sorted.bam"

  samtools view -h "$bam_initial" | python3 -c '
import os, sys, re

# v36 perfectly tiled boundary protection
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

tiled_bounds = {
    "Amp1":  (7729, 7837),
    "Amp2":  (7838, 7940),
    "Amp3":  (7941, 8113),
    "Amp4":  (8114, 8196),
    "Amp5":  (8197, 8293),
    "Amp6":  (8294, 8429),
    "Amp7":  (8430, 8536),
    "Amp8":  (8537, 8638),
    "Amp9":  (8639, 8717),
    "Amp10": (8718, 8877)
}

cigar_re = re.compile(r"(\d+)([MIDNSHPX=])")
# Edge protect is forced to 0 at junctions by using tiled_bounds directly

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
    
    # --- V36: MAXIMUM OVERLAP AMPLICON ASSIGNMENT ---
    best_amp = None
    max_overlap = -1
    
    for name, a_start, a_end in amplicons:
        # Calculate how many bases overlap between the read and this original amplicon
        ovl_start = max(pos, a_start)
        ovl_end = min(ref_end, a_end)
        overlap = ovl_end - ovl_start
        
        if overlap > max_overlap:
            max_overlap = overlap
            best_amp = name
            
    t_start, t_end = tiled_bounds[best_amp]
    # Perfect tiling in 0-based conversion to completely eliminate any gap and overlap
    keep_start = t_start - 1
    keep_end = t_end

    if ref_end <= keep_start or pos >= keep_end:
        continue # Entire read is outside the target amplicon

    ref_p = pos
    read_p = 0
    left_clip = 0
    right_clip = 0
    kept_cigar = []
    new_pos = -1

    for n_str, op in cigar_re.findall(cigar):
        n = int(n_str)
        if op in "M=X":
            start_ref = ref_p
            end_ref = ref_p + n
            ovl_start = max(start_ref, keep_start)
            ovl_end = min(end_ref, keep_end)
            if ovl_start < ovl_end:
                if new_pos == -1: new_pos = ovl_start
                if start_ref < ovl_start: left_clip += (ovl_start - start_ref)
                kept_cigar.append([ovl_end - ovl_start, op])
                if end_ref > ovl_end: right_clip += (end_ref - ovl_end)
            else:
                if end_ref <= keep_start: left_clip += n
                else: right_clip += n
            ref_p += n
            read_p += n
        elif op in "IS":
            if ref_p <= keep_start: left_clip += n
            elif ref_p >= keep_end: right_clip += n
            else: kept_cigar.append([n, op])
            read_p += n
        elif op in "DN":
            start_ref = ref_p
            end_ref = ref_p + n
            ovl_start = max(start_ref, keep_start)
            ovl_end = min(end_ref, keep_end)
            if ovl_start < ovl_end:
                if new_pos == -1: new_pos = ovl_start
                kept_cigar.append([ovl_end - ovl_start, op])
            ref_p += n

    # Strip leading/trailing deletions from the active CIGAR zone
    while kept_cigar and kept_cigar[0][1] in "DN":
        new_pos += kept_cigar[0][0]
        kept_cigar.pop(0)
    while kept_cigar and kept_cigar[-1][1] in "DN":
        kept_cigar.pop()

    if not kept_cigar:
        continue

    # Condense adjacent identical operations
    condensed = []
    for length, op in kept_cigar:
        if not condensed:
            condensed.append([length, op])
        elif condensed[-1][1] == op:
            condensed[-1][0] += length
        else:
            condensed.append([length, op])

    # Reconstruct the Soft-Clipped CIGAR
    final_cigar = ""
    if left_clip > 0: final_cigar += f"{left_clip}S"
    for length, op in condensed:
        final_cigar += f"{length}{op}"
    if right_clip > 0: final_cigar += f"{right_clip}S"

    parts[3] = str(new_pos + 1)
    parts[5] = final_cigar

    try:
        sys.stdout.write("\t".join(parts) + "\n")
    except BrokenPipeError:
        sys.exit(0)
' | samtools view -Sb - | samtools sort -@ "$threads" -o "$bam_trimmed"

  run_log samtools index "$bam_trimmed"

  # Trimmed Variant Calling
  vcf_raw_trimmed="$sample_out/${base}_trimmed_raw.vcf"
  vcf_qual_trimmed="$sample_out/${base}_trimmed_qual_filtered.vcf"
  run_log bash -lc "{ bcftools mpileup -a FORMAT/AD,FORMAT/DP -d $PILEUP_MAX_DEPTH -Q$BASEQ_MIN -q$MAPQ_MIN -Ou -f '$ref' '$bam_trimmed' | bcftools call -mv --ploidy 1 -Ov -o '$vcf_raw_trimmed'; }"
  run_log bcftools filter -i "QUAL>$SAFETY_QUAL && INFO/DP>=$MIN_DEPTH" -Ov -o "$vcf_qual_trimmed" "$vcf_raw_trimmed"

  # Trimmed Annotation
  if [[ -s "$vcf_qual_trimmed" ]]; then
    vcf_annotated_trimmed="$sample_out/${base}_trimmed_annotated_all.vcf"
    run_log bcftools annotate -a "$sample_bed" -c CHROM,FROM,TO,RegionType -h "$sample_out/hdr.txt" -Ov -o "$vcf_annotated_trimmed" "$vcf_qual_trimmed"

    run_log bash -lc "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_trimmed' | bcftools filter -e 'POS=7898 || POS=7899 || POS=8595' -Ov -o '$sample_out/${base}_trimmed_snps.vcf'"
    run_log bash -lc "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_trimmed' | bcftools filter -e 'RegionType=\"HP_Region\" || RegionType=\"Blacklist_Site\"' -Ov -o '$sample_out/${base}_trimmed_clean.vcf'"
    run_log bcftools filter -i 'RegionType="HP_Region"' -Ov -o "$sample_out/${base}_trimmed_homopolymers.vcf" "$vcf_annotated_trimmed"
    
    rm -f "$sample_bed" "$sample_out/hdr.txt"
  fi

  # =========================================================================
  # DIAGNOSTIC DEPTH COMPARISON
  # =========================================================================
  echo "[Diagnostic] Calculating average depth across genome..." | tee -a "$log_file"
  depth_baseline=$(samtools depth -a -q "$BASEQ_MIN" "$bam_baseline" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
  depth_trimmed=$(samtools depth -a -q "$BASEQ_MIN" "$bam_trimmed" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
  
  echo "---------------------------------------------------" | tee -a "$log_file"
  echo "DIAGNOSTIC SUMMARY: $base" | tee -a "$log_file"
  echo "Baseline Average Depth : $depth_baseline" | tee -a "$log_file"
  echo "Trimmed Average Depth  : $depth_trimmed" | tee -a "$log_file"
  echo "---------------------------------------------------" | tee -a "$log_file"

  run_log fastqc "$bam_trimmed" -o "$sample_out"
  echo "Completed Diagnostic Split: $base" | tee -a "$log_file"
done

echo "All outputs written to: $run_out"
