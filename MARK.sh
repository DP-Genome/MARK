#!/usr/bin/env bash
# MARK.sh
# tag: MARK
#
# PURPOSE
# ---------------------------------------------------------------------------
# Diagnostic Split-Track Pipeline for Nanopore mitochondrial data.
# ---------------------------------------------------------------------------

set -euo pipefail

# Data files (reference, regions BED, adapter lists) ship in the same folder as this
# script. Find that folder, following symlinks, so a run started from any directory
# uses the pipeline's own files rather than whatever happens to sit in the cwd.
_self="${BASH_SOURCE[0]}"
while [[ -L "$_self" ]]; do
  _dir="$(cd -P "$(dirname "$_self")" && pwd)"
  _self="$(readlink "$_self")"
  if [[ "$_self" != /* ]]; then _self="$_dir/$_self"; fi
done
SCRIPT_DIR="$(cd -P "$(dirname "$_self")" && pwd)"
unset _self _dir
# A bare file name that ships beside this script means that copy. Paths, and names the
# script folder does not have, are used exactly as given.
mark_resolve() {
  if [[ -n "$1" && "$1" != */* && -f "$SCRIPT_DIR/$1" ]]; then printf '%s\n' "$SCRIPT_DIR/$1"
  else printf '%s\n' "$1"; fi
}

VERSION="1.2.5"

if [[ $# -lt 1 || "$1" == "-h" || "$1" == "--help" ]]; then
  echo -e "MARK Pipeline (ONT) v$VERSION"
  echo -e "Usage: MARK.sh <input_fastq_file_or_folder>\n"
  echo -e "DESCRIPTION:"
  echo -e "  Diagnostic Split-Track Pipeline for Nanopore mitochondrial data."
  echo -e "  This script relies on environment variables for configuration."
  echo -e "\nCOMMON OVERRIDE VARIABLES (with defaults):"
  echo -e "  threads=\"8\"              Number of CPU threads to use"
  echo -e "  MIN_DEPTH=\"10\"           Minimum depth for variant calling"
  echo -e "  QS_MIN=\"10\"              Minimum average read quality (fastplong)"
  echo -e "  CUTADAPT_STALL_SECS=\"600\" Seconds without cutadapt output before it is killed and rerun single-core"
  echo -e "  MIN_LEN=\"90\"             Minimum read length before trimming"
  echo -e "  MAX_LEN=\"1500\"           Maximum read length (LENSAFE filter)"
  echo -e "  MIN_LEN_POST=\"90\"        Minimum read length after trimming"
  echo -e "  MAX_LEN_POST=\"300\"       Maximum read length after trimming"
  echo -e "  EXTRA_TRIM=\"0\"           Extra bases to trim from both ends"
  echo -e "  ref=\"linearized_mtdna.fasta\"     Reference FASTA file (default: the copy beside this script)"
  echo -e "  regions_bed=\"linearized_regions.bed\" Amplicon BED file (default: the copy beside this script)"
  echo -e "  RUN_NAME=\"\"             Name for the output folder (default: auto-generated)"
  echo -e "  OUTPUT_DIR=\"\"           Where to create it (default: beside the input)"
  echo -e "\nEXAMPLE EXECUTIONS:"
  echo -e "  # Run with defaults:"
  echo -e "  MARK.sh /path/to/fastqs"
  echo -e "\n  # Run overriding threads and minimum depth:"
  echo -e "  threads=16 MIN_DEPTH=20 MARK.sh /path/to/fastqs"
  exit 1
fi

input_path="$1"
[[ -e "$input_path" ]] || { echo "Error: Input path '$input_path' not found."; exit 1; }

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
CUTADAPT_STALL_SECS="${CUTADAPT_STALL_SECS:-600}"
CUTADAPT_ERR="${CUTADAPT_ERR:-0.10}" 
CUTADAPT_OVL="${CUTADAPT_OVL:-5}"    

ref="${ref:-linearized_mtdna.fasta}"
ref="$(mark_resolve "$ref")"
mmi_index="${mmi_index:-${ref}.mmi}"
regions_bed="${regions_bed:-linearized_regions.bed}"
regions_bed="$(mark_resolve "$regions_bed")"

SAFETY_QUAL="${SAFETY_QUAL:-20}"
STRICT_QUAL="${STRICT_QUAL:-60}"
MIN_DEPTH="${MIN_DEPTH:-10}"       
PILEUP_MAX_DEPTH="${PILEUP_MAX_DEPTH:-100000}"
# bcftools keeps a SEPARATE depth cap for indel candidates (--max-idepth,
# default 250). Raising -d alone leaves it at 250, which silently stops indel
# candidate generation above 250x - i.e. everywhere in this panel. Tie it to
# PILEUP_MAX_DEPTH so one setting governs both and nothing is capped below it.
# The snps/clean outputs stay SNP-only (they filter on TYPE="snp"); this only
# restores indels to the annotated_all / qual_filtered review VCFs, which is
# where they are meant to appear.
PILEUP_MAX_IDEPTH="${PILEUP_MAX_IDEPTH:-$PILEUP_MAX_DEPTH}"
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
  main_dir="$(dirname "$input_path")"
  if [[ -f "$SCRIPT_DIR/MARK_Adapter_List_ONT.txt" ]]; then
    ADAPTER_FILE="$SCRIPT_DIR/MARK_Adapter_List_ONT.txt"
  elif [[ -f "$main_dir/MARK_Adapter_List_ONT.txt" ]]; then
    ADAPTER_FILE="$main_dir/MARK_Adapter_List_ONT.txt"
  elif [[ -f "MARK_Adapter_List_ONT.txt" ]]; then
    ADAPTER_FILE="$(pwd)/MARK_Adapter_List_ONT.txt"
  fi
fi
ADAPTER_FILE="$(mark_resolve "$ADAPTER_FILE")"

# --- Exact amplicon trimming (v1.2.4) ------------------------------------------------
# TRIM_MODE=adapters : validated path - trim with the adapter lists.
# TRIM_MODE=primers  : trim every read at its own amplicon primers. A read is kept only if it
#                      carries the forward AND reverse primer of the SAME amplicon, with nothing
#                      but library adapter outside them, and is cut exactly at the outer primer
#                      ends. The primers stay on the read on purpose: the product coordinates in
#                      CRM_Nested_products.bed include them, and the primer-footprint step later
#                      soft-clips them in the BAM, so no primer base reaches the caller.
TRIM_MODE="${TRIM_MODE:-adapters}"
[[ "$TRIM_MODE" == "adapters" || "$TRIM_MODE" == "primers" ]] || { echo "Error: TRIM_MODE must be adapters or primers (got '$TRIM_MODE')."; exit 1; }
PRIMERS_BED="${PRIMERS_BED:-CRM_Nested_primers_empirical.bed}"
PRIMERS_BED="$(mark_resolve "$PRIMERS_BED")"
PRIMER_ERR="${PRIMER_ERR:-0.15}"
PRIMER_MIN_OVERLAP="${PRIMER_MIN_OVERLAP:-18}"
# A low primer-pair share never stops the run (v1.2.5). Reads without a primer pair are always
# removed; a library below PRIMER_WARN_PCT is flagged in WARNINGS.txt, its log and the run summary,
# so a kit or primer mismatch cannot pass unnoticed.
PRIMER_WARN_PCT="${PRIMER_WARN_PCT:-10}"
# Rejected reads kept for review in <library>_primer_rejects_review*.fastq: up to this many per step.
PRIMER_DIAG_READS="${PRIMER_DIAG_READS:-20000}"
if [[ "$TRIM_MODE" == "primers" ]]; then
  [[ -f "$PRIMERS_BED" ]] || { echo "Error: primer BED '$PRIMERS_BED' not found (needed for TRIM_MODE=primers)."; exit 1; }
fi

# Write the cutadapt primer files for this run from the BED and the reference, so they can never
# drift from the coordinates the rest of the pipeline uses.
build_primer_files() {
  python3 - "$PRIMERS_BED" "$ref" "$1" <<'PYEOF'
import sys
bed_path, ref_path, out = sys.argv[1:4]
seq = "".join(l.strip() for l in open(ref_path) if not l.startswith(">"))
bed = {}
for line in open(bed_path):
    if line.startswith("#") or not line.strip():
        continue
    f = line.split("\t")
    bed[f[3]] = (int(f[1]), int(f[2]))          # BED: 0-based, half open
amps = ["Amp%d" % i for i in range(1, 11)]
rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
FWD = {a: seq[bed[a + "_FWD"][0]:bed[a + "_FWD"][1]] for a in amps}
REV = {a: seq[bed[a + "_REV"][0]:bed[a + "_REV"][1]] for a in amps}   # + strand: end of the product
linked = open(out + "/primers_linked.fa", "w")
anch = open(out + "/primers_anchored.fa", "w")
r1 = open(out + "/primers_read1.fa", "w")
r2 = open(out + "/primers_read2.fa", "w")
single = open(out + "/primers_single.fa", "w")
for a in amps:
    linked.write(">%s\n%s...%s\n" % (a, FWD[a], REV[a]))
    anch.write(">%s\n^%s...%s$\n" % (a, FWD[a], REV[a]))
    # --pair-adapters: entry n of the read 1 list pairs with entry n of the read 2 list
    r1.write(">%s_F\n^%s\n>%s_R\n^%s\n" % (a, FWD[a], a, rc(REV[a])))
    r2.write(">%s_R\n^%s\n>%s_F\n^%s\n" % (a, rc(REV[a]), a, FWD[a]))
    single.write(">%s_F\n%s\n>%s_R\n%s\n" % (a, FWD[a], a, rc(REV[a])))
for fh in (linked, anch, r1, r2, single):
    fh.close()
PYEOF
}

# reads kept per amplicon, from the amplicon name cutadapt appended to each read name
primer_amplicon_counts() {   # $1 = fastq, $2 = label, $3 = output tsv
  awk -v label="$2" 'NR%4==1 { n=split($1, p, "_"); a=p[n]; sub(/_[FR]$/, "", a);
                       if (a !~ /^Amp[0-9]+$/) a="unassigned"; c[a]++ }
       END { for (k in c) printf "%s\t%s\t%d\n", label, k, c[k] }' "$1" | sort >> "$3"
}

# Why were reads rejected? Classified on at most PRIMER_DIAG_READS of them, so the run summary can
# separate "no primer at all", "only one primer" and "primers from two different amplicons".
# The classified reads are also appended, each labelled with its reason, to the review fastq ($5).
primer_reject_reasons() {    # $1 = rejected fastq, $2 = primer dir, $3 = out tsv, $4 = label, $5 = review fastq
  local n; n=$(count_fastq_reads "$1")
  if [[ "$n" -eq 0 ]]; then printf "%s\treject_total\t0\n" "$4" >> "$3"; return; fi
  head -n $((PRIMER_DIAG_READS * 4)) "$1" > "$1.diag.fq"
  cutadapt -g "file:$2/primers_single.fa" --action=none --revcomp -e "$PRIMER_ERR" --times 2 \
    --info-file "$1.diag.info" -o /dev/null "$1.diag.fq" > /dev/null 2>&1 || true
  python3 - "$1.diag.info" "$n" "$4" "$1.diag.fq" "$5" <<'PYEOF' >> "$3"
import collections, sys
info, total, label, diag_fq, review = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
hits, seen = {}, []
for line in open(info):
    f = line.rstrip("\n").split("\t")
    if len(f) < 8:              # no primer matched: cutadapt writes a short line (name, -1, seq, qual)
        name = f[0].split()[0]
        if name not in hits:
            seen.append(name)
            hits[name] = set()
        continue
    name, err, adapter = f[0].split()[0], f[-10], f[-4]
    if name not in hits:
        seen.append(name)
        hits[name] = set()
    if err != "-1" and adapter != "no_adapter":
        hits[name].add(adapter.split("_")[0])
c, why = collections.Counter(), {}
for name in seen:
    amps = hits[name]
    if not amps:
        kind, detail = "no primer found", ""
    elif len(amps) > 1:
        kind, detail = "primers from two amplicons", " (%s)" % "+".join(sorted(amps))
    else:
        kind, detail = "only one primer of one amplicon", " (%s)" % next(iter(amps))
    c[kind] += 1
    why[name] = kind + detail
with open(diag_fq) as fh, open(review, "a") as out:
    while True:
        h = fh.readline()
        if not h:
            break
        rest = fh.readline() + fh.readline() + fh.readline()
        name = h[1:].split()[0]
        out.write("@%s mark_reject=%s: %s\n" % (name, label, why.get(name, "not classified")) + rest)
scale = total / max(len(seen), 1)
print("%s\treject_total\t%d" % (label, total))
for reason, k in c.most_common():
    print("%s\treject: %s\t%d" % (label, reason, round(k * scale)))
PYEOF
  rm -f "$1.diag.fq" "$1.diag.info"
}

# A library whose primer-pair share is below PRIMER_WARN_PCT: never stop, but make it impossible to
# miss - a banner on screen and in the log, and a line in WARNINGS.txt at the top of the run folder.
flag_low_primer_share() {    # $1 = library, $2 = pct, $3 = kept, $4 = total, $5 = review file(s)
  local w="$run_out/WARNINGS.txt"
  [[ -f "$w" ]] || printf "# MARK v%s warnings for this run. Libraries listed here were processed, but check them.\nlibrary\tprimer_pair_pct\twarn_below_pct\tkept\tof\treview_reads\n" "$VERSION" > "$w"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$PRIMER_WARN_PCT" "$3" "$4" "$5" >> "$w"
  {
    echo "################################################################################"
    echo "## WARNING: $1 - only $2% of reads carry a primer pair of one amplicon"
    echo "##          (flag level ${PRIMER_WARN_PCT}%). The library WAS processed with the $3 reads that do."
    echo "##          Check the kit and primer BED; rejected reads for review: $5"
    echo "################################################################################"
  } | tee -a "$log_file" >&2
}

# Nanopore reads still carry the sequencing adapter, so the primer match cannot be anchored.
# Keep a read only if what lay OUTSIDE its two primers was library adapter (or nothing): this is
# what rejects read-through products, which carry primers of two different amplicons.
ont_flank_filter() {   # $1 = reads into cutadapt, $2 = primer-trimmed reads, $3 = output, $4 = stats tsv, $5 = review fastq
  python3 - "$1" "$2" "$3" "$4" "$5" "$PRIMER_DIAG_READS" <<'PYEOF'
import sys
src, trimmed, out, stats, review, cap = sys.argv[1:7]
cap = int(cap)
L_JUNC, R_JUNC = "CTTCCGATCT", "AGATCGGAAG"     # adapter bases next to the amplicon on each side
MAX_EDITS, MIN_FLANK = 3, 6

def edits_at_end(flank, motif):                 # best alignment of motif to the END of flank
    t = flank[-(len(motif) + 3):]
    n, m = len(motif), len(t)
    prev = [0] * (m + 1)            # a prefix of the flank may be skipped for free
    for i in range(1, n + 1):
        cur = [i] + [0] * m
        for j in range(1, m + 1):
            cur[j] = min(prev[j - 1] + (motif[i - 1] != t[j - 1]), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return prev[m]                  # the motif must reach the end of the flank

def edits_at_start(flank, motif):
    return edits_at_end(flank[:len(motif) + 3][::-1], motif[::-1])

rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def reads(path):
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            yield h, fh.readline().rstrip("\n"), fh.readline(), fh.readline()

# A retained product must not contain library adapter INSIDE it: that is a chimera whose
# junction happens to sit between two primer matches.
INTERNAL = ["AATGATACGGCGACCACCGAGATCTACAC", "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
            "AGATCGGAAGAGC", "ATCTCGTATGCCGTCTTCTGCTTG"]
INTERNAL = INTERNAL + [x[::-1].translate(str.maketrans("ACGTN", "TGCAN")) for x in INTERNAL]

kept = rejected = chimeric = 0
src_iter = reads(src)
src_name, src_seq, src_qual = None, None, None
rev = open(review, "a")
n_rev = 0
def to_review(reason):
    global n_rev
    if n_rev < cap:
        rev.write("@%s mark_reject=adapter_check: %s\n%s\n+\n%s" % (src_name, reason, src_seq, src_qual))
        n_rev += 1
with open(out, "w") as fh:
    for h, seq, plus, qual in reads(trimmed):
        name = h[1:].split()[0]
        base = name.rsplit("_Amp", 1)[0]
        while src_name != base:
            try:
                sh, sseq, _, squal = next(src_iter)
            except StopIteration:
                break
            src_name, src_seq, src_qual = sh[1:].split()[0], sseq, squal
        if src_name != base:
            continue
        hay = src_seq
        i = hay.find(seq)
        if i < 0:
            hay = rc(src_seq)
            i = hay.find(seq)
        if i < 0:                                # cutadapt corrected bases; nothing to check
            fh.write(h + seq + "\n" + plus + qual); kept += 1; continue
        left, right = hay[:i], hay[i + len(seq):]
        lok = len(left) < MIN_FLANK or edits_at_end(left, L_JUNC) <= MAX_EDITS
        rok = len(right) < MIN_FLANK or edits_at_start(right, R_JUNC) <= MAX_EDITS
        if not (lok and rok):
            rejected += 1; to_review("non-adapter sequence outside the primers")
        elif any(m in seq for m in INTERNAL):
            chimeric += 1; to_review("adapter inside the product (chimera)")
        else:
            fh.write(h + seq + "\n" + plus + qual); kept += 1
rev.close()
with open(stats, "a") as fh:
    fh.write(f"flank_filter\tkept\t{kept}\n"
             f"flank_filter\treject: non-adapter sequence outside the primers\t{rejected}\n"
             f"flank_filter\treject: adapter inside the product (chimera)\t{chimeric}\n")
PYEOF
}

MINIMAP_EXTRA="${MINIMAP_EXTRA:-}"
if [[ "$TRIM_MODE" == "primers" && -z "$MINIMAP_EXTRA" ]]; then MINIMAP_EXTRA="--end-bonus 15"; fi

if [[ -z "$ADAPTER_FILE" || ! -f "$ADAPTER_FILE" ]]; then
  echo "Error: adapter file not found"
  exit 1
fi

command -v fastplong >/dev/null 2>&1 || { echo "Error: fastplong not found."; exit 1; }

# --- Output Folder Setup ---
run_ts="$(date +%Y%m%d_%H%M%S)"
input_name="$(basename "$input_path")"
default_run_name="${pipeline_name}_${input_name}_${run_ts}_output"

# OUTPUT_DIR is the directory the run folder is created in. It defaults to the
# directory holding the input, so a run sits BESIDE its input rather than inside
# it, and does so regardless of where the script was invoked from.
if [[ -n "${OUTPUT_DIR:-}" ]]; then
  out_base_dir="$OUTPUT_DIR"
  mkdir -p "$out_base_dir" 2>/dev/null || true
elif [[ -d "$input_path" ]]; then
  out_base_dir="$(dirname "$(cd "$input_path" && pwd)")"
else
  out_base_dir="$(cd "$(dirname "$input_path")" && pwd)"
fi
if [[ ! -d "$out_base_dir" ]]; then
  echo "Error: output directory '$out_base_dir' does not exist and could not be created."
  exit 1
fi
out_base_dir="$(cd "$out_base_dir" && pwd)"

# RUN_NAME lets the launcher (or the user) name the output folder explicitly.
# It is reduced to a single, safe path component; anything unusable falls back
# to the default naming convention.
if [[ -n "${RUN_NAME:-}" ]]; then
  run_name="$(basename "$RUN_NAME")"
  run_name="${run_name//[^A-Za-z0-9._-]/_}"
  run_name="${run_name#.}"
  [[ -n "$run_name" ]] || run_name="$default_run_name"
else
  run_name="$default_run_name"
fi
run_out="$out_base_dir/$run_name"

# Never write a new run into a folder that already holds results.
if [[ -d "$run_out" && -n "$(ls -A "$run_out" 2>/dev/null)" ]]; then
  echo "Error: output folder '$run_out' already exists and is not empty."
  echo "       Pick a different run name, or move the existing folder aside."
  exit 1
fi
mkdir -p "$run_out"

ref_chrom_name=$(head -n1 "$ref" | cut -d ' ' -f1 | tr -d '>')

run_log() {
  {
    printf 'Running: %s\n' "$*"
    "$@"
    printf '\n'
  } 2>&1 | tee -a "$log_file"
}

# cutadapt's multi-core mode can deadlock on macOS: the parent process waits forever
# on worker processes that have already died, nothing times out, and the run hangs
# with no error. So watch the output file instead: if it stops growing for
# CUTADAPT_STALL_SECS, kill the run and repeat it single-core, which has no worker
# processes to lose. The output is identical either way, because cutadapt keeps reads
# in input order in multi-core mode. A second stall is a hard error, never a skip.
_cutadapt_out_size() { stat -c %s "$1" 2>/dev/null || stat -f %z "$1" 2>/dev/null || echo 0; }
run_cutadapt() {
  local args=("$@") out="" i attempt pid rc size last idle stalled
  for ((i = 0; i < ${#args[@]}; i++)); do
    if [[ "${args[$i]}" == "-o" ]]; then out="${args[$((i + 1))]}"; fi
  done
  for attempt in 1 2; do
    if [[ -n "$out" ]]; then rm -f "$out"; fi
    printf 'Running: cutadapt %s\n' "${args[*]}" >> "$log_file"
    cutadapt "${args[@]}" >> "$log_file" 2>&1 &
    pid=$!
    last=-1; idle=0; stalled=0
    while kill -0 "$pid" 2>/dev/null; do
      sleep 1
      size=$(_cutadapt_out_size "$out")
      if [[ "$size" != "$last" ]]; then last=$size; idle=0; else idle=$((idle + 1)); fi
      if (( idle >= CUTADAPT_STALL_SECS )) && kill -0 "$pid" 2>/dev/null; then
        stalled=1
        pkill -9 -P "$pid" 2>/dev/null || true
        kill -9 "$pid" 2>/dev/null || true
        break
      fi
    done
    if wait "$pid" 2>/dev/null; then rc=0; else rc=$?; fi
    if (( stalled == 0 )); then
      printf '\n' >> "$log_file"
      return "$rc"
    fi
    echo "[cutadapt] output stopped growing for ${CUTADAPT_STALL_SECS}s; killed and rerunning single-core" | tee -a "$log_file" >&2
    for ((i = 0; i < ${#args[@]}; i++)); do
      if [[ "${args[$i]}" == "--cores" ]]; then args[$((i + 1))]=1; fi
    done
  done
  echo "Error: cutadapt stalled again, single-core, writing $out" | tee -a "$log_file" >&2
  return 1
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
  echo "## Pipeline Version: $VERSION"
  echo "## Input Path     : $(realpath "$input_path")"
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
  echo "## PILEUP_IDEPTH  : $PILEUP_MAX_IDEPTH"
  echo "##"
  echo "## --- VARIANT FILTERING ---"
  echo "## MIN_DEPTH      : $MIN_DEPTH"
  echo "## SAFETY_QUAL    : $SAFETY_QUAL"
  echo "## STRICT_QUAL    : $STRICT_QUAL"
  echo "##"
  echo "## =============================================================================="
  printf "sample\tstage\treads_in\treads_out\tdropped_reads\tdropped_pct\n"
} > "$run_summary_file"

if [[ -d "$input_path" ]]; then
  files=("$input_path"/*.fastq "$input_path"/*.fastq.gz)
else
  files=("$input_path")
fi

primer_dir="$run_out/primer_files"
if [[ "$TRIM_MODE" == "primers" ]]; then
  mkdir -p "$primer_dir"
  build_primer_files "$primer_dir" || { echo "Error: could not build the primer files from $PRIMERS_BED"; exit 1; }
  echo "[Primers] built cutadapt primer files from $(basename "$PRIMERS_BED") in $primer_dir"
fi

for fq in "${files[@]}"; do
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
  
  # Adapter removal, both ends, one pass.
  #
  # This replaces a two-pass "-a then -g". These libraries are Illumina
  # amplicons read end-to-end on ONT, so the adapter sits at BOTH ends of the
  # molecule. "-a" treats every match as a 3-prime adapter and deletes the match
  # and everything after it, so a match at the 5-prime end deleted the entire
  # amplicon, and --minimum-length then discarded the read.
  #
  # "-b" treats a match that includes the first base of the read as a 5-prime
  # adapter (removing it and keeping what follows) and any other match as a
  # 3-prime adapter. --times 2 lets one read lose an adapter from each end.
  t5="$sample_out/${base}_qsGE${QS_MIN}_adaptertrim.fastq"
  if [[ "$TRIM_MODE" == "primers" ]]; then
    # Cut at the amplicon's own primers. --action=retain keeps the two primer matches and removes
    # everything outside them, which is exactly the adapter, index and flow-cell sequence.
    amp_stats="$sample_out/${base}_amplicon_trim.tsv"
    : > "$amp_stats"
    primer_raw="$sample_out/${base}_primertrim_raw.fastq"
    primer_rej="$sample_out/${base}_primertrim_rejected.fastq"
    run_cutadapt -g "file:$primer_dir/primers_linked.fa" --revcomp --action=retain \
      --rename '{id}_{adapter_name}' --error-rate "$PRIMER_ERR" --overlap "$PRIMER_MIN_OVERLAP" \
      --untrimmed-output "$primer_rej" --cores "$threads" -o "$primer_raw" "$qs_fq" >/dev/null
    primer_raw_reads=$(count_fastq_reads "$primer_raw")
    printf "%s\t2_PrimerTrim\t%s\t%s\t%s\t%s%%\n" "$base" "$qs_reads" "$primer_raw_reads" \
      "$((qs_reads - primer_raw_reads))" \
      "$(awk -v d=$((qs_reads - primer_raw_reads)) -v i="$qs_reads" 'BEGIN{ if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')" \
      >> "$run_summary_file"
    review="$sample_out/${base}_primer_rejects_review.fastq"
    : > "$review"
    ont_flank_filter "$qs_fq" "$primer_raw" "$t5" "$amp_stats" "$review"
    primer_reject_reasons "$primer_rej" "$primer_dir" "$amp_stats" "primer_step" "$review"
    primer_amplicon_counts "$t5" "kept" "$amp_stats"
    t5_reads=$(count_fastq_reads "$t5")
    printf "%s\t2b_FlankFilter\t%s\t%s\t%s\t%s%%\n" "$base" "$primer_raw_reads" "$t5_reads" \
      "$((primer_raw_reads - t5_reads))" \
      "$(awk -v d=$((primer_raw_reads - t5_reads)) -v i="$primer_raw_reads" 'BEGIN{ if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')" \
      >> "$run_summary_file"
    pct=$(awk -v k="$t5_reads" -v i="$qs_reads" 'BEGIN{ if(i>0) printf "%.1f", (k/i)*100; else print "0.0" }')
    echo "[Primers] $t5_reads of $qs_reads reads ($pct%) carry both primers of one amplicon" | tee -a "$log_file"
    share_step="2c_PrimerPairShare"
    if awk -v p="$pct" -v w="$PRIMER_WARN_PCT" 'BEGIN{ exit (p+0 < w+0) ? 0 : 1 }'; then
      share_step="2c_PrimerPairShare_LOW_see_WARNINGS"
      flag_low_primer_share "$base" "$pct" "$t5_reads" "$qs_reads" "$(basename "$review")"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s%%\n" "$base" "$share_step" "$qs_reads" "$t5_reads" "$((qs_reads - t5_reads))" \
      "$(awk -v k="$t5_reads" -v i="$qs_reads" 'BEGIN{ if(i>0) printf "%.2f", (1-k/i)*100; else print "0.00" }')" >> "$run_summary_file"
    rm -f "$primer_raw" "$primer_rej"
  else
    run_cutadapt -b "file:$ADAPTER_FILE" --times 2 --error-rate "$CUTADAPT_ERR" --overlap "$CUTADAPT_OVL" --minimum-length "$MIN_LEN" --cores "$threads" -o "$t5" "$qs_fq" >/dev/null
    t5_reads=$(count_fastq_reads "$t5")
  fi
  t5_drop=$((qs_reads - t5_reads))
  t5_pct=$(awk -v d="$t5_drop" -v i="$qs_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t2_AdapterTrim\t%s\t%s\t%s\t%s%%\n" "$base" "$qs_reads" "$t5_reads" "$t5_drop" "$t5_pct" >> "$run_summary_file"
  
  final_fq="$sample_out/${base}_qsGE${QS_MIN}_trim5_u${EXTRA_TRIM}x2.fastq"
  if [[ "$EXTRA_TRIM" -gt 0 ]]; then
    run_cutadapt -u "$EXTRA_TRIM" -u "-$EXTRA_TRIM" --minimum-length "$MIN_LEN" --cores "$threads" -o "$final_fq" "$t5" >/dev/null
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
  # With exact primer trimming the read ends are real template, so let the aligner reach them
  # instead of soft-clipping an error-containing primer end (--end-bonus).
  run_log bash -c "{ minimap2 -ax map-ont $MINIMAP_EXTRA -t $threads '$mmi_index' '$final_fq' | samtools view -Sb - | samtools sort -@ $threads -o '$bam_initial'; }"
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
  run_log bash -c "{ bcftools mpileup -a FORMAT/AD,FORMAT/DP -d $PILEUP_MAX_DEPTH --max-idepth $PILEUP_MAX_IDEPTH -Q$BASEQ_MIN -q$MAPQ_MIN -Ou -f '$ref' '$bam_baseline' | bcftools call -mv --ploidy 1 -Ov -o '$vcf_raw_base'; }"
  run_log bcftools filter -i "QUAL>$SAFETY_QUAL && INFO/DP>=$MIN_DEPTH" -Ov -o "$vcf_qual_base" "$vcf_raw_base"

  # Baseline Annotation
  if [[ -s "$vcf_qual_base" ]]; then
    sample_bed="$sample_out/${base}_targets.bed"
    sed "s/PLACEHOLDER_CHROM/$ref_chrom_name/g" "$regions_bed" > "$sample_bed"
    echo '##INFO=<ID=RegionType,Number=1,Type=String,Description="Type of region (HP_Region or Blacklist_Site)">' > "$sample_out/hdr.txt"

    vcf_annotated_base="$sample_out/${base}_baseline_annotated_all.vcf"
    run_log bcftools annotate -a "$sample_bed" -c CHROM,FROM,TO,RegionType -h "$sample_out/hdr.txt" -Ov -o "$vcf_annotated_base" "$vcf_qual_base"

    run_log bash -c "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_base' | bcftools filter -e 'POS=7898 || POS=7899 || POS=8595' -Ov -o '$sample_out/${base}_baseline_snps.vcf'"
    run_log bash -c "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_base' | bcftools filter -e 'RegionType=\"HP_Region\" || RegionType=\"Blacklist_Site\"' -Ov -o '$sample_out/${base}_baseline_clean.vcf'"
    run_log bcftools filter -i 'RegionType="HP_Region"' -Ov -o "$sample_out/${base}_baseline_homopolymers.vcf" "$vcf_annotated_base"
  fi


  # =========================================================================
  # TRACK B: EXPERIMENTAL (PYTHON TRIMMING/SOFT-CLIPPING)
  # =========================================================================
  echo "[Pipeline] Generating TRIMMED track..." | tee -a "$log_file"
  bam_trimmed="$sample_out/${base}_trimmed_sorted.bam"

  samtools view -h "$bam_initial" | python3 -c '
import os, sys, re

# Perfectly tiled boundary protection
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

# Callable insert per amplicon: the product minus its own primer footprints.
# These are the ten inserts published for the PowerSeq CRM Nested System
# (Vinueza-Espinosa et al. 2023, Electrophoresis 44:1423-1434). Unlike the midpoint
# tiles they replace, adjacent inserts overlap - a base covered by two amplicons is
# sequenced by two independent molecules, so both are counted. Each fragment is
# still counted once. Assignment stays on these same spans; assigning on the full
# product span instead lets a neighbour with a long primer arm capture reads it
# then clips away, losing real coverage.
tiled_bounds = {
    "Amp1":  (7729, 7842),
    "Amp2":  (7832, 7941),
    "Amp3":  (7939, 8124),
    "Amp4":  (8103, 8202),
    "Amp5":  (8190, 8315),
    "Amp6":  (8271, 8437),
    "Amp7":  (8421, 8542),
    "Amp8":  (8531, 8649),
    "Amp9":  (8627, 8721),
    "Amp10": (8714, 8877)
}

# Product spans and their primer footprints. Promega does not publish the primers,
# so these come from the fixed read termini seen when sequencing this kit
# (see CRM_Nested_primers_empirical.bed).
#   name, product_start, product_end, forward_primer_end, reverse_primer_start  (1-based)
products = [
    ("Amp1",  7702, 7868, 7728, 7843),
    ("Amp2",  7810, 7964, 7831, 7942),
    ("Amp3",  7913, 8149, 7938, 8125),
    ("Amp4",  8079, 8225, 8102, 8203),
    ("Amp5",  8166, 8337, 8189, 8316),
    ("Amp6",  8258, 8465, 8270, 8438),
    ("Amp7",  8394, 8578, 8420, 8543),
    ("Amp8",  8503, 8674, 8530, 8650),
    ("Amp9",  8608, 8745, 8626, 8722),
    ("Amp10", 8687, 8904, 8713, 8878),
]
PRIMER_TOL = 10


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
    
    # --- PRIMER-FOOTPRINT TRIMMING ---
    # A read is no longer assigned to one amplicon and cut to that amplicon bounds.
    # Instead each end of the read is checked against the known product termini:
    # if this read starts where a product starts, its own forward primer is removed;
    # if it ends where a product ends, its own reverse primer is removed. Everything
    # between is template and is kept, whichever amplicons it spans.
    #
    # This matters for shared and hybrid fragments. A fragment running from the Amp4
    # forward primer to the Amp5 reverse primer physically covers both inserts; the
    # max-overlap rule gave it entirely to Amp5 and clipped away its Amp4 half,
    # discarding real template.
    #
    # An end that matches no product terminus is left alone - that is a read whose
    # end is genuine template (an unmerged mate, or a partial read), not primer.
    start1 = pos + 1
    end1 = ref_end
    keep_start = pos
    keep_end = ref_end
    matched = False
    for name, p_start, p_end, f_end, r_start in products:
        if abs(start1 - p_start) <= PRIMER_TOL:
            matched = True
            keep_start = max(keep_start, f_end)
            # a read of this product whose far end stops inside the products own
            # reverse primer is showing primer sequence there, not template
            if r_start <= end1 <= p_end + PRIMER_TOL:
                keep_end = min(keep_end, r_start - 1)
        if abs(end1 - p_end) <= PRIMER_TOL:
            matched = True
            keep_end = min(keep_end, r_start - 1)
            if p_start - PRIMER_TOL <= start1 <= f_end:
                keep_start = max(keep_start, f_end)

    if not matched:
        # Neither end lines up with a known product, so we cannot tell which primers
        # this fragment carries. Fall back to the conservative primer-aware rule:
        # assign to the amplicon it overlaps most and clip to that insert. Such reads
        # are rare but tend to occur where genuine coverage is thin, so leaving them
        # untrimmed would let primer sequence stand in for missing template.
        best_amp = None
        max_overlap = -1
        for name, a_start, a_end in amplicons:
            overlap = min(ref_end, a_end) - max(pos, a_start)
            if overlap > max_overlap:
                max_overlap = overlap
                best_amp = name
        if best_amp is None or max_overlap <= 0:
            continue
        t_start, t_end = tiled_bounds[best_amp]
        keep_start = t_start - 1
        keep_end = t_end

    # off-target guard: the read must still touch at least one amplicon insert
    if not any(min(keep_end, a_end) - max(keep_start, a_start) > 0
               for _, a_start, a_end in amplicons):
        continue

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
  run_log bash -c "{ bcftools mpileup -a FORMAT/AD,FORMAT/DP -d $PILEUP_MAX_DEPTH --max-idepth $PILEUP_MAX_IDEPTH -Q$BASEQ_MIN -q$MAPQ_MIN -Ou -f '$ref' '$bam_trimmed' | bcftools call -mv --ploidy 1 -Ov -o '$vcf_raw_trimmed'; }"
  run_log bcftools filter -i "QUAL>$SAFETY_QUAL && INFO/DP>=$MIN_DEPTH" -Ov -o "$vcf_qual_trimmed" "$vcf_raw_trimmed"

  # Trimmed Annotation
  if [[ -s "$vcf_qual_trimmed" ]]; then
    vcf_annotated_trimmed="$sample_out/${base}_trimmed_annotated_all.vcf"
    run_log bcftools annotate -a "$sample_bed" -c CHROM,FROM,TO,RegionType -h "$sample_out/hdr.txt" -Ov -o "$vcf_annotated_trimmed" "$vcf_qual_trimmed"

    run_log bash -c "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_trimmed' | bcftools filter -e 'POS=7898 || POS=7899 || POS=8595' -Ov -o '$sample_out/${base}_trimmed_snps.vcf'"
    run_log bash -c "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_trimmed' | bcftools filter -e 'RegionType=\"HP_Region\" || RegionType=\"Blacklist_Site\"' -Ov -o '$sample_out/${base}_trimmed_clean.vcf'"
    run_log bcftools filter -i 'RegionType="HP_Region"' -Ov -o "$sample_out/${base}_trimmed_homopolymers.vcf" "$vcf_annotated_trimmed"
    
# =========================================================================
# AUTOMATED HETEROPLASMY EXTRACTION (VAF)
# =========================================================================
echo "[Processing] Extracting heteroplasmy summary for $base..." | tee -a "$log_file"

# Define the output path strictly inside the sample's timestamped folder
het_output="$sample_out/${base}_final_sample_report.tsv"

python3 -c '
import sys
import os
import glob

sample_dir = sys.argv[1]
base_name = sys.argv[2]
het_output = os.path.join(sample_dir, f"{base_name}_final_sample_report.tsv")

def get_region(pos):
    if 16024 <= pos <= 16365:
        return "Control Region (HVI)", "16024-16365"
    elif 73 <= pos <= 340:
        return "Control Region (HVII)", "73-340"
    elif 438 <= pos <= 574:
        return "Control Region (HVIII)", "438-574"
    elif (15990 <= pos <= 16569) or (1 <= pos <= 700):
        return "Control Region", "Control Region"
    else:
        return "Coding Region", "Whole Genome"

main_keywords = ["snps", "clean", "annotated_all", "indels", "qual_filtered"]

main_variants = {}
other_variants = {}

vcf_files = glob.glob(os.path.join(sample_dir, "*.vcf"))

def parse_vcf(vcf_path):
    variants = {}
    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            
            chrom = parts[0]
            lin_pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            fmt = parts[8].split(":")
            sample_data = parts[9].split(":")
            
            try:
                ad_idx = fmt.index("AD")
                dp_idx = fmt.index("DP")
                
                ad_field = sample_data[ad_idx]
                tot_dp = int(sample_data[dp_idx]) if sample_data[dp_idx] != "." else 0
                
                if tot_dp == 0:
                    continue
                    
                ad_parts = ad_field.split(",")
                ref_dp = int(ad_parts[0]) if ad_parts[0] != "." else 0
                alt_dp = int(ad_parts[1]) if len(ad_parts) > 1 and ad_parts[1] != "." else 0
                
                vaf = (alt_dp / tot_dp) * 100.0 if tot_dp > 0 else 0.0
                
                if 1.0 <= vaf <= 100.0:
                    variant_key = f"{chrom}_{lin_pos}_{ref}_{alt}"
                    variants[variant_key] = {
                        "chrom": chrom,
                        "lin_pos": lin_pos,
                        "ref": ref,
                        "alt": alt,
                        "vaf": vaf
                    }
            except (ValueError, IndexError, ZeroDivisionError):
                continue
    return variants

for vcf_file in vcf_files:
    fname = os.path.basename(vcf_file).lower()
    is_main = any(kw in fname for kw in main_keywords)
    
    parsed_vars = parse_vcf(vcf_file)
    for k, v in parsed_vars.items():
        if is_main:
            main_variants[k] = v
        else:
            other_variants[k] = v

# Remove duplicates from other_variants that are already in main_variants
for k in main_variants:
    if k in other_variants:
        del other_variants[k]

with open(het_output, "w") as out:
    def write_variants(var_dict):
        sorted_vars = sorted(var_dict.values(), key=lambda x: x["lin_pos"])
        for v in sorted_vars:
            lin_pos = v["lin_pos"]
            ref = v["ref"]
            alt = v["alt"]
            vaf = v["vaf"]
            
            real_pos = (lin_pos + 8284) % 16569
            if real_pos == 0: real_pos = 16569
            
            region_name, range_cov = get_region(real_pos)
            
            disp_ref = ref
            variant_called = ""
            notes = ""
            het_status = ""
            
            if len(ref) == 1 and len(alt) == 1:
                disp_ref = ref
                variant_called = f"m.{real_pos}{alt}"
                notes = f"m.{real_pos}{ref}>{alt}"
                het_status = f"Homoplasmic {vaf:.0f}%" if vaf > 95.0 else f"{vaf:.0f}% PHP"
            elif len(ref) > len(alt):
                del_seq = ref[len(alt):]
                disp_ref = del_seq
                del_start = real_pos + len(alt)
                if len(del_seq) == 1:
                    variant_called = f"m.{del_start}del"
                else:
                    del_end = del_start + len(del_seq) - 1
                    variant_called = f"m.{del_start}_{del_end}del"
                notes = "Deletion"
                het_status = f"Homoplasmic {vaf:.0f}%" if vaf > 95.0 else f"{vaf:.0f}% LHP"
            else:
                disp_ref = "-"
                ins_seq = alt[len(ref):]
                variant_called = f"m.{real_pos}.1{ins_seq}"
                notes = "Insertion"
                het_status = f"Homoplasmic {vaf:.0f}%" if vaf > 95.0 else f"{vaf:.0f}% LHP"
                
            out.write(f"{region_name}\t{range_cov}\t{disp_ref}\t{variant_called}\t{het_status}\t{notes}\n")

    out.write("### MAIN RESULTS ###\n")
    out.write("Region / Locus\tRange Covered\trCRS Reference Nucleotide\tVariant Called\tHeteroplasmy % (VAF)\tNotes / HGVS Equivalent\n")
    write_variants(main_variants)
    
    out.write("\n### OTHER OBSERVATIONS ###\n")
    out.write("Region / Locus\tRange Covered\trCRS Reference Nucleotide\tVariant Called\tHeteroplasmy % (VAF)\tNotes / HGVS Equivalent\n")
    write_variants(other_variants)

' "$sample_out" "$base"

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

if [[ -f "$run_out/WARNINGS.txt" ]]; then
  {
    echo "################################################################################"
    echo "## THIS RUN HAS WARNINGS - see $run_out/WARNINGS.txt"
    awk -F'\t' 'NR>2 { printf "##   %s: %s%% of reads carry a primer pair\n", $1, $2 }' "$run_out/WARNINGS.txt"
    echo "################################################################################"
  } >&2
fi
echo "All outputs written to: $run_out"
