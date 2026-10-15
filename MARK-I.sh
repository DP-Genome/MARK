#!/usr/bin/env bash
# MARK-I.sh
# tag: MARK-I
#
# PURPOSE
# ---------------------------------------------------------------------------
# Illumina Merged-Read Diagnostic Split-Track Pipeline for mitochondrial data.
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

VERSION="1.4.1"

if [[ $# -lt 1 || "$1" == "-h" || "$1" == "--help" ]]; then
  echo -e "MARK Pipeline (Illumina) v$VERSION"
  echo -e "Usage: MARK-I.sh <input_fastq_file_or_folder>\n"
  echo -e "DESCRIPTION:"
  echo -e "  Illumina Merged-Read Diagnostic Split-Track Pipeline for mitochondrial data."
  echo -e "  This script relies on environment variables for configuration."
  echo -e "\nCOMMON OVERRIDE VARIABLES (with defaults):"
  echo -e "  threads=\"8\"              Number of CPU threads to use"
  echo -e "  MIN_DEPTH=\"10\"           Minimum depth for variant calling"
  echo -e "  READQ=\"20\"               Quality score requirement (fastp)"
  echo -e "  CUTADAPT_STALL_SECS=\"600\" Seconds without cutadapt output before it is killed and rerun single-core"
  echo -e "  MIN_LEN=\"90\"             Min RAW fragment length (validated floor; dropped)"
  echo -e "  MIN_LEN_POST=\"90\"        Shortest product kept after trimming (same window as nanopore)"
  echo -e "  MAX_LEN_POST=\"300\"       Longest product kept after trimming (same window as nanopore)"
  echo -e "  MIN_ALIGNED_POST=\"30\"    Template bases a read must align after primer clipping (else read and mate dropped)"
  echo -e "  ALIGNED_WARN_PCT=\"20\"    Flag an amplicon in WARNINGS.txt if that check removes this % of its fragments"
  echo -e "  (MIN_DEPTH also flags an amplicon in WARNINGS.txt when any insert position is below it in the final BAM)"
  echo -e "  EXTRA_TRIM=\"0\"           Extra bases to trim from both ends"
  echo -e "  ref=\"linearized_mtdna.fasta\"     Reference FASTA file (default: the copy beside this script)"
  echo -e "  regions_bed=\"linearized_regions.bed\" Amplicon BED file (default: the copy beside this script)"
  echo -e "  RUN_NAME=\"\"             Name for the output folder (default: auto-generated)"
  echo -e "  OUTPUT_DIR=\"\"           Where to create it (default: beside the input)"
  echo -e "\nEXAMPLE EXECUTIONS:"
  echo -e "  # Run with defaults:"
  echo -e "  MARK-I.sh /path/to/fastqs"
  echo -e "\n  # Run overriding threads, minimum quality, and max read length:"
  echo -e "  threads=16 READQ=25 MIN_DEPTH=20 MARK-I.sh /path/to/fastqs"
  exit 1
fi

input_path="$1"
[[ -e "$input_path" ]] || { echo "Error: Input path '$input_path' not found."; exit 1; }

# --- Configuration ---
threads="${threads:-8}"
pipeline_name="MARK-I"

DISCARD_WARN_PCT="${DISCARD_WARN_PCT:-5}"

# fastp settings
READQ="${READQ:-20}"
UNQUAL_PCT="${UNQUAL_PCT:-40}"
MIN_LEN="${MIN_LEN:-90}"
MIN_LEN_POST="${MIN_LEN_POST:-90}"
MAX_LEN_POST="${MAX_LEN_POST:-300}"
# v1.2.7: after primer clipping, a read needs this many template bases aligned (M/=/X), or it is
# dropped together with its mate. Catches reads that pass every length filter but align only a
# few bases, e.g. an Illumina read 2 that turns to poly-G right after its primer (27S3M121S).
MIN_ALIGNED_POST="${MIN_ALIGNED_POST:-30}"
# v1.2.8: an amplicon that loses this share of its fragments (%) to that check, and at least
# ALIGNED_WARN_MIN of them, is flagged in WARNINGS.txt. The library is still processed.
ALIGNED_WARN_PCT="${ALIGNED_WARN_PCT:-20}"
ALIGNED_WARN_MIN="${ALIGNED_WARN_MIN:-10}"
N_BASE_LIMIT="${N_BASE_LIMIT:-5}"

# cutadapt settings
CUTADAPT_ERR="${CUTADAPT_ERR:-0.10}"
CUTADAPT_OVL="${CUTADAPT_OVL:-5}"
EXTRA_TRIM="${EXTRA_TRIM:-0}"
CUTADAPT_STALL_SECS="${CUTADAPT_STALL_SECS:-600}"

ref="${ref:-linearized_mtdna.fasta}"
ref="$(mark_resolve "$ref")"
regions_bed="${regions_bed:-linearized_regions.bed}"
regions_bed="$(mark_resolve "$regions_bed")"

# HIERARCHICAL FILTERS
SAFETY_QUAL="${SAFETY_QUAL:-20}"
STRICT_QUAL="${STRICT_QUAL:-60}"
MIN_DEPTH="${MIN_DEPTH:-10}"       
# v1.2.10: no limit by default (0 = no limit; was 100,000). A pileup that reaches a limit skips reads and
# drops real calls (see flag_pileup_cap). bcftools gets its own maximum for 0, never 0 itself.
PILEUP_MAX_DEPTH="${PILEUP_MAX_DEPTH:-0}"
BCF_MAX_DEPTH="$PILEUP_MAX_DEPTH"; [[ "$BCF_MAX_DEPTH" -le 0 ]] && BCF_MAX_DEPTH=2147483647
# v1.4: variants are called with freebayes (local realignment, one copy). PILEUP_MAX_DEPTH > 0 is passed
# to it as --limit-coverage; 0 (default) = no limit.
_fb_user="${FREEBAYES:-}"                 # set by the user: used as given, never replaced
FREEBAYES="${FREEBAYES:-freebayes}"
# freebayes is installed in its own conda environment (its bioconda builds cannot share one with the
# validated htslib/fastp versions on Linux and Apple Silicon). If it is not on PATH and FREEBAYES is not set,
# look for an environment named mark_fb beside the active one or in the usual conda locations.
if [[ -z "$_fb_user" ]] && ! command -v "$FREEBAYES" >/dev/null 2>&1; then
  _ce="${CONDA_EXE:-}"
  for _fb in "${CONDA_PREFIX:-/nonexistent}/../mark_fb/bin/freebayes" "${CONDA_PREFIX:-/nonexistent}/envs/mark_fb/bin/freebayes" \
             "${_ce%/bin/conda}/envs/mark_fb/bin/freebayes" "$HOME/miniforge3/envs/mark_fb/bin/freebayes" \
             "$HOME/miniconda3/envs/mark_fb/bin/freebayes" "$HOME/anaconda3/envs/mark_fb/bin/freebayes" \
             "/opt/homebrew/anaconda3/envs/mark_fb/bin/freebayes" "/opt/anaconda3/envs/mark_fb/bin/freebayes"; do
    if [[ -x "$_fb" ]]; then FREEBAYES="$_fb"; break; fi
  done
  unset _ce _fb
fi
unset _fb_user
FB_LIMIT=""; [[ "$PILEUP_MAX_DEPTH" -gt 0 ]] && FB_LIMIT="--limit-coverage $PILEUP_MAX_DEPTH"
# bcftools keeps a SEPARATE depth cap for indel candidates (--max-idepth,
# default 250). Raising -d alone leaves it at 250, which silently stops indel
# candidate generation above 250x - i.e. everywhere in this panel. Tie it to
# PILEUP_MAX_DEPTH so one setting governs both and nothing is capped below it.
# The snps/clean outputs stay SNP-only (they filter on TYPE="snp"); this only
# restores indels to the annotated_all / qual_filtered review VCFs, which is
# where they are meant to appear.
PILEUP_MAX_IDEPTH="${PILEUP_MAX_IDEPTH:-$BCF_MAX_DEPTH}"; [[ "$PILEUP_MAX_IDEPTH" -le 0 ]] && PILEUP_MAX_IDEPTH=2147483647
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
  if [[ -f "$SCRIPT_DIR/MARK_Adapter_List_Illumina.txt" ]]; then
    ADAPTER_FILE="$SCRIPT_DIR/MARK_Adapter_List_Illumina.txt"
  elif [[ -f "$main_dir/MARK_Adapter_List_Illumina.txt" ]]; then
    ADAPTER_FILE="$main_dir/MARK_Adapter_List_Illumina.txt"
  elif [[ -f "MARK_Adapter_List_Illumina.txt" ]]; then
    ADAPTER_FILE="$(pwd)/MARK_Adapter_List_Illumina.txt"
  fi
fi
ADAPTER_FILE="$(mark_resolve "$ADAPTER_FILE")"

# --- Trimming mode --------------------------------------------------------------------
# Both modes remove the library adapters, and in both the primer-footprint step later soft-clips
# the primer bases in the BAM. They differ in how a read is cut and which reads are kept:
# TRIM_MODE=standard : the pre-v1.2.4 path - adapter sequences from the adapter list are cut
#                      wherever they are found, and any read long enough afterwards is kept.
# TRIM_MODE=exact    : the default since v1.2.6 - every read is cut at its own amplicon
#                      primers. A read is kept only if it carries the forward AND reverse primer of
#                      the SAME amplicon, with nothing but library adapter outside them, and is cut
#                      exactly at the outer primer ends. The primers stay on the read on purpose: the
#                      product coordinates in CRM_Nested_products.bed include them, and the
#                      primer-footprint step soft-clips them in the BAM, so no primer base reaches
#                      the caller.
# The v1.2.4-v1.2.6 names still work: adapters = standard, primers = exact.
TRIM_MODE="${TRIM_MODE:-exact}"
case "$TRIM_MODE" in
  adapters) echo "[Note] TRIM_MODE=adapters is now called TRIM_MODE=standard"; TRIM_MODE=standard ;;
  primers)  echo "[Note] TRIM_MODE=primers is now called TRIM_MODE=exact"; TRIM_MODE=exact ;;
esac
[[ "$TRIM_MODE" == "standard" || "$TRIM_MODE" == "exact" ]] || { echo "Error: TRIM_MODE must be standard or exact (got '$TRIM_MODE')."; exit 1; }
PRIMERS_BED="${PRIMERS_BED:-CRM_Nested_primers_empirical.bed}"
PRIMERS_BED="$(mark_resolve "$PRIMERS_BED")"
# v1.2.9: the amplicon inserts, for the per-amplicon depth check (flagged in WARNINGS.txt)
INSERTS_BED="${INSERTS_BED:-CRM_Nested_inserts.bed}"
INSERTS_BED="$(mark_resolve "$INSERTS_BED")"
# v1.2.9: the primer-clipping step reads both BEDs, in either mode
[[ -f "$PRIMERS_BED" && -f "$INSERTS_BED" ]] || { echo "Error: primer BED '$PRIMERS_BED' or inserts BED '$INSERTS_BED' not found."; exit 1; }
PRIMER_ERR="${PRIMER_ERR:-0.10}"
PRIMER_MIN_OVERLAP="${PRIMER_MIN_OVERLAP:-3}"
# A low primer-pair share never stops the run (v1.2.5). Reads without a primer pair are always
# removed; a library below PRIMER_WARN_PCT is flagged in WARNINGS.txt, its log and the run summary,
# so a kit or primer mismatch cannot pass unnoticed.
PRIMER_WARN_PCT="${PRIMER_WARN_PCT:-25}"
# Rejected reads kept for review in <library>_primer_rejects_review*.fastq: up to this many per step.
PRIMER_DIAG_READS="${PRIMER_DIAG_READS:-20000}"
if [[ "$TRIM_MODE" == "exact" ]]; then
  [[ -f "$PRIMERS_BED" ]] || { echo "Error: primer BED '$PRIMERS_BED' not found (needed for TRIM_MODE=exact)."; exit 1; }
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
# One line per warning in WARNINGS.txt at the top of the run folder (created on the first warning).
mark_warn() {   # $1 library, $2 check, $3 amplicon, $4 value %, $5 flag level %, $6 count, $7 of, $8 detail
  local w="$run_out/WARNINGS.txt"
  [[ -f "$w" ]] || printf "# MARK v%s warnings for this run. Libraries listed here were processed, but check them.\nlibrary\tcheck\tamplicon\tvalue\tflag_level\tcount\tof\tdetail\n" "$VERSION" > "$w"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$@" >> "$w"
}

flag_low_primer_share() {    # $1 = library, $2 = pct, $3 = kept, $4 = total, $5 = review file(s)
  mark_warn "$1" "primer-pair share below flag level" "all" "$2%" "$PRIMER_WARN_PCT%" "$3" "$4" "fragments kept; rejected reads for review: $5"
  {
    echo "################################################################################"
    echo "## WARNING: $1 - only $2% of fragments carry a primer pair of one amplicon"
    echo "##          (flag level ${PRIMER_WARN_PCT}%). The library WAS processed with the $3 fragments that do."
    echo "##          Check the kit and primer BED; rejected reads for review: $5"
    echo "################################################################################"
  } | tee -a "$log_file" >&2
}

# v1.2.9: an amplicon with any insert position below MIN_DEPTH in the final BAM. Depth is counted as the
# caller counts it, except that a read with a deletion at a position covers it (samtools depth -J), so a
# homozygous deletion is not taken for a gap. Unclassified barcode bins are not samples and are skipped.
flag_low_depth() {   # $1 = library, $2 = final BAM
  [[ -f "$INSERTS_BED" ]] || return 0
  local tsv="$sample_out/${1}_amplicon_depth.tsv"
  samtools depth -a -J -q "$BASEQ_MIN" -Q "$MAPQ_MIN" -b "$INSERTS_BED" "$2" 2>/dev/null | \
    awk -F'\t' -v bed="$INSERTS_BED" -v floor="$MIN_DEPTH" '
      BEGIN { OFS = "\t"
              while ((getline l < bed) > 0) { if (l ~ /^#/ || l == "") continue
                split(l, f, "\t"); n++; s[n] = f[2]; e[n] = f[3]; a[n] = f[4]; sub(/_insert$/, "", a[n]); mn[n] = -1 } }
      { for (i = 1; i <= n; i++) if ($2 > s[i] && $2 <= e[i]) {
          if (mn[i] < 0 || $3 < mn[i]) { mn[i] = $3; at[i] = $2 }
          tot[i] += $3; c[i]++; if ($3 < floor) lo[i]++ } }
      END { print "amplicon", "positions", "mean_depth", "min_depth", "min_at_rcrs", "positions_below_min_depth"
            for (i = 1; i <= n; i++) print a[i], c[i] + 0, (c[i] ? sprintf("%.0f", tot[i] / c[i]) : 0), (mn[i] < 0 ? 0 : mn[i]),
                                          (at[i] ? (at[i] + 8283) % 16569 + 1 : "-"), lo[i] + 0 }' > "$tsv"
  [[ "$1" == *unclassified* ]] && return 0
  awk -F'\t' 'NR > 1 && ($6 > 0 || $2 == 0)' "$tsv" | while IFS=$'\t' read -r amp npos mean mind at below; do
    [[ "$npos" -eq 0 ]] && { below="all"; mind=0; at="-"; }
    mark_warn "$1" "amplicon depth below MIN_DEPTH" "$amp" "${mind}x" "${MIN_DEPTH}x" "$below" "$npos" \
      "insert positions below ${MIN_DEPTH}x in the final BAM, lowest at rCRS $at; see $(basename "$tsv")"
    {
      echo "################################################################################"
      echo "## WARNING: $1 $amp - $below of its $npos insert positions are below ${MIN_DEPTH}x in the final BAM"
      echo "##          (lowest ${mind}x at rCRS $at). Variants there cannot be called."
      echo "##          Per-amplicon depth: $(basename "$tsv")"
      echo "################################################################################"
    } | tee -a "$log_file" >&2
  done
}

# v1.2.10: the variant caller reads at most PILEUP_MAX_DEPTH reads per position (v1.4: freebayes
# --limit-coverage). When a limit is reached, reads are dropped and real calls can be lost.
# The highest depth in the final BAM (deletions counted, no quality filter, as bcftools counts reads) must stay
# below it; otherwise the library is flagged and PILEUP_MAX_DEPTH must be raised.
flag_pileup_cap() {   # $1 = library, $2 = final BAM
  local maxd
  maxd=$(samtools depth -J "$2" 2>/dev/null | awk 'BEGIN { m = 0 } $3 > m { m = $3; p = $2 } END { printf "%d\t%s", m, (p ? (p + 8283) % 16569 + 1 : "-") }')
  local mdepth="${maxd%%$'\t'*}" mat="${maxd##*$'\t'}"
  echo "[Depth] highest read depth in the final BAM: ${mdepth}x at rCRS $mat (PILEUP_MAX_DEPTH ${PILEUP_MAX_DEPTH}; 0 = no limit)" | tee -a "$log_file"
  if [[ "$PILEUP_MAX_DEPTH" -gt 0 && "${mdepth:-0}" -ge "$PILEUP_MAX_DEPTH" ]]; then
    mark_warn "$1" "depth at or above PILEUP_MAX_DEPTH" "all" "${mdepth}x" "${PILEUP_MAX_DEPTH}x" "-" "-" \
      "the variant caller drops reads beyond PILEUP_MAX_DEPTH, so calls may be missing; rerun with PILEUP_MAX_DEPTH above ${mdepth}"
    {
      echo "################################################################################"
      echo "## WARNING: $1 - read depth reaches ${mdepth}x (rCRS $mat), at or above PILEUP_MAX_DEPTH (${PILEUP_MAX_DEPTH})."
      echo "##          The variant caller drops reads beyond it, so variant calls may be MISSING."
      echo "##          Rerun with PILEUP_MAX_DEPTH set above ${mdepth}."
      echo "################################################################################"
    } | tee -a "$log_file" >&2
  fi
}

# v1.2.8: an amplicon that lost ALIGNED_WARN_PCT % or more of its fragments to the post-alignment check.
flag_aligned_losses() {   # $1 = library, $2 = per-amplicon table
  [[ -s "$2" ]] || return 0
  awk -F'\t' -v p="$ALIGNED_WARN_PCT" -v m="$ALIGNED_WARN_MIN" 'NR>1 && $3+0 >= m+0 && $4+0 >= p+0' "$2" |
  while IFS=$'\t' read -r amp total removed pct; do
    mark_warn "$1" "fragments removed after alignment above flag level" "$amp" "$pct%" "$ALIGNED_WARN_PCT%" "$removed" "$total" \
      "fragments with a read that aligns fewer than $MIN_ALIGNED_POST template bases, or no mate; see $(basename "$2")"
    {
      echo "################################################################################"
      echo "## WARNING: $1 $amp - $removed of $total fragments ($pct%) were removed after alignment:"
      echo "##          one of their reads aligns fewer than $MIN_ALIGNED_POST template bases, or has no aligned mate"
      echo "##          (flag level ${ALIGNED_WARN_PCT}%). Coverage of $amp comes from the fragments left."
      echo "##          Check the review BAM; per-amplicon counts: $(basename "$2")"
      echo "################################################################################"
    } | tee -a "$log_file" >&2
  done
}

# Second merge pass: pairs that fastp could not merge, quality-trimmed from the 3' end first.
# (--cut_right in the FIRST pass is worse on every amplicon, so it is only used here.)
MERGE2_CUT_Q="${MERGE2_CUT_Q:-20}"
MERGE2_CUT_WINDOW="${MERGE2_CUT_WINDOW:-4}"
MERGE2_MIN_LEN="${MERGE2_MIN_LEN:-30}"
BWA_EXTRA="${BWA_EXTRA:-}"

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

command -v "$FREEBAYES" >/dev/null 2>&1 || { echo "Error: freebayes (the variant caller since v1.4) was not found. Install it in its own environment:"; echo "  conda create -n mark_fb -c conda-forge -c bioconda freebayes=1.3.10"; echo "MARK finds an environment named mark_fb by itself; otherwise set FREEBAYES=/path/to/freebayes."; exit 1; }
command -v fastp >/dev/null 2>&1 || { echo "Error: fastp not found."; exit 1; }

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
  echo "## TRIM_MODE      : $TRIM_MODE$([[ "$TRIM_MODE" == "exact" ]] && echo " (reads cut to one whole amplicon)" || echo " (adapter-list trimming)")"
  echo "## FASTP MERGE    : READQ=$READQ | UNQUAL_PCT=$UNQUAL_PCT | MIN_LEN(raw)=$MIN_LEN"
  echo "## MIXED MATCH    : ERR=$CUTADAPT_ERR | OVL=$CUTADAPT_OVL | MIN_LEN_POST=$MIN_LEN_POST "
  echo "## EXTRA_TRIM     : $EXTRA_TRIM"
  echo "## MAX_LEN_POST   : $MAX_LEN_POST"
  echo "## MIN_LEN_POST   : $MIN_LEN_POST"
  echo "##"
  echo "## --- ALIGNMENT & CALLING ---"
  echo "## MAPQ_MIN       : $MAPQ_MIN"
  echo "## BASEQ_MIN      : $BASEQ_MIN"
  echo "## PILEUP_DEPTH   : $PILEUP_MAX_DEPTH"
  echo "## VARIANT CALLER : freebayes $("$FREEBAYES" --version 2>/dev/null | awk '{print $NF}') ($(command -v "$FREEBAYES"))"
  echo "## PILEUP_IDEPTH  : $PILEUP_MAX_IDEPTH"
  echo "## MIN_ALIGNED_POST: $MIN_ALIGNED_POST (template bases a read must align after primer clipping; else read and mate removed)"
  echo "## ALIGNED_WARN   : ${ALIGNED_WARN_PCT}% of an amplicon's fragments removed by that check, at least $ALIGNED_WARN_MIN (flagged in WARNINGS.txt)"
  echo "## DEPTH_WARN     : an amplicon with any insert position below MIN_DEPTH ($MIN_DEPTH x) in the final BAM (flagged in WARNINGS.txt)"
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
  files=(
    "$input_path"/*_R1*.fastq "$input_path"/*_R1*.fastq.gz
    "$input_path"/*.R1*.fastq "$input_path"/*.R1*.fastq.gz
  )
else
  files=("$input_path")
fi

primer_dir="$run_out/primer_files"
if [[ "$TRIM_MODE" == "exact" ]]; then
  mkdir -p "$primer_dir"
  build_primer_files "$primer_dir" || { echo "Error: could not build the primer files from $PRIMERS_BED"; exit 1; }
  echo "[Primers] built cutadapt primer files from $(basename "$PRIMERS_BED") in $primer_dir"
fi

for r1 in "${files[@]}"; do
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
  input_reads=$(count_fastq_reads "$r1") 

  # v1.2.9: no length-safety pass on Illumina. NextSeq reads are at most 151 bases, so a maximum length can
  # never apply; empty reads (nothing left after the sequencer's adapter trimming) fail the 90-base floor
  # at merging instead. Nanopore keeps its MAX_LEN pass.
  lensafe_r1="$r1"; lensafe_r2="$r2"
  len_reads=$input_reads

  # --- STEP 1: Merge & Filter ---
  merged_fq="$sample_out/${base}_merged.fastq"
  # Only merged reads are aligned. Pairs fastp cannot merge are written to their own
  # files so they can be counted and inspected, but they are not used.
  unmerged_r1="$sample_out/${base}_unmerged_R1.fastq"
  unmerged_r2="$sample_out/${base}_unmerged_R2.fastq"
  # v1.2.9: a fragment is kept only when read 1 and read 2 both pass and merge into one fragment.
  # When one mate fails filtering, fastp discards the pair, and no mate is rescued (either mode).
  fp_html="$sample_out/${base}_fastp_merge.html"
  fp_json="$sample_out/${base}_fastp_merge.json"

  run_log fastp \
    -i "$lensafe_r1" -I "$lensafe_r2" \
    --merge --trim_poly_g \
    --merged_out "$merged_fq" \
    --out1 "$unmerged_r1" --out2 "$unmerged_r2" \
    --qualified_quality_phred "$READQ" \
    --unqualified_percent_limit "$UNQUAL_PCT" \
    --length_required "$MIN_LEN" \
    --n_base_limit "$N_BASE_LIMIT" \
    --thread "$threads" \
    --html "$fp_html" \
    --json "$fp_json" 2> /dev/null

  merged_se=$(count_fastq_reads "$merged_fq")
  unmerged_pairs=$(count_fastq_reads "$unmerged_r1")
  merged_reads=$merged_se
  if [[ "$TRIM_MODE" == "exact" ]]; then
    echo "[Merge] $unmerged_pairs pairs did not merge in the first pass; they get a second pass" | tee -a "$log_file"
  else
    echo "[Merge] $unmerged_pairs pairs could not be merged and are not used" | tee -a "$log_file"
  fi
  merged_drop=$((len_reads * 2 - merged_reads))
  merged_pct=$(awk -v d="$merged_drop" -v i="$((len_reads * 2))" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t1_MergeFilter(Reads)\t%s\t%s\t%s\t%s%%\n" "$base" "$((len_reads * 2))" "$merged_reads" "$merged_drop" "$merged_pct" >> "$run_summary_file"

  # --- STEP 2: Preprocessing ---
  run_log fastqc "$merged_fq" -o "$sample_out" > /dev/null 2>&1

  if [[ "$TRIM_MODE" == "exact" ]]; then
    amp_stats="$sample_out/${base}_amplicon_trim.tsv"
    : > "$amp_stats"
    # --- pass 2: merge the leftovers after 3' quality trimming ---
    merged2_fq="$sample_out/${base}_merged_pass2.fastq"
    still_r1="$sample_out/${base}_unmerged2_R1.fastq"
    still_r2="$sample_out/${base}_unmerged2_R2.fastq"
    if [[ -s "$unmerged_r1" ]]; then
      run_log fastp -i "$unmerged_r1" -I "$unmerged_r2" \
        --merge --merged_out "$merged2_fq" --out1 "$still_r1" --out2 "$still_r2" \
        --cut_right --cut_right_window_size "$MERGE2_CUT_WINDOW" --cut_right_mean_quality "$MERGE2_CUT_Q" \
        --disable_quality_filtering --length_required "$MERGE2_MIN_LEN" --thread "$threads" \
        --html "$sample_out/${base}_fastp_merge2.html" --json "$sample_out/${base}_fastp_merge2.json" 2> /dev/null
    else
      : > "$merged2_fq"; : > "$still_r1"; : > "$still_r2"
    fi
    merged2_n=$(count_fastq_reads "$merged2_fq")
    still_n=$(count_fastq_reads "$still_r1")
    printf "%s\t1b_MergePass2\t%s\t%s\t%s\t%s%%\n" "$base" "$unmerged_pairs" "$merged2_n" \
      "$((unmerged_pairs - merged2_n))" \
      "$(awk -v d=$((unmerged_pairs - merged2_n)) -v i="$unmerged_pairs" 'BEGIN{ if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')" \
      >> "$run_summary_file"

    # --- merged reads of both passes: must BE one whole product, primer to primer ---
    merged_all="$sample_out/${base}_merged_all.fastq"
    cat "$merged_fq" "$merged2_fq" > "$merged_all"
    merged_all_n=$(count_fastq_reads "$merged_all")
    merged_exact="$sample_out/${base}_merged_exact.fastq"
    merged_rej="$sample_out/${base}_merged_rejected.fastq"
    run_cutadapt -g "file:$primer_dir/primers_anchored.fa" --revcomp --action=retain \
      --rename '{id}_{adapter_name}' --error-rate "$PRIMER_ERR" \
      --untrimmed-output "$merged_rej" --cores "$threads" -o "$merged_exact" "$merged_all" >/dev/null
    merged_exact_n=$(count_fastq_reads "$merged_exact")
    printf "%s\t2_PrimerTrim(merged)\t%s\t%s\t%s\t%s%%\n" "$base" "$merged_all_n" "$merged_exact_n" \
      "$((merged_all_n - merged_exact_n))" \
      "$(awk -v d=$((merged_all_n - merged_exact_n)) -v i="$merged_all_n" 'BEGIN{ if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')" \
      >> "$run_summary_file"
    review_m="$sample_out/${base}_primer_rejects_review_merged.fastq"
    : > "$review_m"
    primer_reject_reasons "$merged_rej" "$primer_dir" "$amp_stats" "merged" "$review_m"
    primer_amplicon_counts "$merged_exact" "kept_merged" "$amp_stats"

    # --- pairs that never merged, in either pass: not used (v1.2.9) ---
    # Only fragments whose two reads merge are kept, so every read in the BAM is one whole amplicon
    # read from both ends. Unmerged pairs are counted here and dropped.
    pe_r1="$sample_out/${base}_pairs_exact_R1.fastq"; pe_r2="$sample_out/${base}_pairs_exact_R2.fastq"
    : > "$pe_r1"; : > "$pe_r2"
    pe_pairs=0
    printf "%s\t2b_UnmergedPairs_dropped\t%s\t0\t%s\t%s%%\n" "$base" "$still_n" "$still_n" \
      "$(awk -v i="$still_n" 'BEGIN{ if(i>0) print "100.00"; else print "0.00" }')" >> "$run_summary_file"
    echo "[Merge] $still_n pairs did not merge in either pass and are not used" | tee -a "$log_file"

    kept_frag=$((merged_exact_n + pe_pairs))
    in_frag=$((merged_all_n + still_n))
    pct=$(awk -v k="$kept_frag" -v i="$in_frag" 'BEGIN{ if(i>0) printf "%.1f", (k/i)*100; else print "0.0" }')
    echo "[Primers] $kept_frag of $in_frag fragments ($pct%) are a merged, complete single-amplicon product" | tee -a "$log_file"
    share_step="2c_PrimerPairShare"
    if awk -v p="$pct" -v w="$PRIMER_WARN_PCT" 'BEGIN{ exit (p+0 < w+0) ? 0 : 1 }'; then
      share_step="2c_PrimerPairShare_LOW_see_WARNINGS"
      flag_low_primer_share "$base" "$pct" "$kept_frag" "$in_frag" "$(basename "$review_m")"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s%%\n" "$base" "$share_step" "$in_frag" "$kept_frag" "$((in_frag - kept_frag))" \
      "$(awk -v k="$kept_frag" -v i="$in_frag" 'BEGIN{ if(i>0) printf "%.2f", (1-k/i)*100; else print "0.00" }')" >> "$run_summary_file"
    # v1.2.9: the same product-length window as nanopore (MIN_LEN_POST-MAX_LEN_POST). A plain length filter,
    # so nothing else about the reads changes. Every whole product is 144-237 bp, so it removes nothing here.
    merged_len="$sample_out/${base}_merged_exact_len.fastq"
    awk -v lo="$MIN_LEN_POST" -v hi="$MAX_LEN_POST" \
      'NR % 4 == 1 { h = $0; next } NR % 4 == 2 { s = $0; next } NR % 4 == 3 { p = $0; next }
       { if (length(s) >= lo + 0 && length(s) <= hi + 0) print h "\n" s "\n" p "\n" $0 }' "$merged_exact" > "$merged_len"
    len_ok=$(count_fastq_reads "$merged_len")
    printf "%s\t2d_ProductLength(%s-%sbp)\t%s\t%s\t%s\t%s%%\n" "$base" "$MIN_LEN_POST" "$MAX_LEN_POST" "$merged_exact_n" "$len_ok" \
      "$((merged_exact_n - len_ok))" \
      "$(awk -v d=$((merged_exact_n - len_ok)) -v i="$merged_exact_n" 'BEGIN{ if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')" \
      >> "$run_summary_file"
    mv "$merged_len" "$merged_exact"
    final_fq="$merged_exact"
    rm -f "$merged_all" "$merged_rej" "$merged2_fq" "$still_r1" "$still_r2"
  else
  
  t3="$sample_out/${base}_trim3.fastq"
  run_cutadapt -a "file:$ADAPTER_FILE" --error-rate "$CUTADAPT_ERR" --overlap "$CUTADAPT_OVL" --minimum-length "$MIN_LEN_POST" --cores "$threads" -o "$t3" "$merged_fq" >/dev/null
  t3_reads=$(count_fastq_reads "$t3")
  t3_drop=$((merged_se - t3_reads))
  t3_pct=$(awk -v d="$t3_drop" -v i="$merged_se" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t2_Trim3(merged)\t%s\t%s\t%s\t%s%%\n" "$base" "$merged_se" "$t3_reads" "$t3_drop" "$t3_pct" >> "$run_summary_file"

  t5="$sample_out/${base}_trim5.fastq"
  run_cutadapt -g "file:$ADAPTER_FILE" --error-rate "$CUTADAPT_ERR" --overlap "$CUTADAPT_OVL" --minimum-length "$MIN_LEN_POST" --cores "$threads" -o "$t5" "$t3" >/dev/null
  t5_reads=$(count_fastq_reads "$t5")
  t5_drop=$((t3_reads - t5_reads))
  t5_pct=$(awk -v d="$t5_drop" -v i="$t3_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t3_Trim5\t%s\t%s\t%s\t%s%%\n" "$base" "$t3_reads" "$t5_reads" "$t5_drop" "$t5_pct" >> "$run_summary_file"
  
  final_fq="$sample_out/${base}_trim5_u${EXTRA_TRIM}x2.fastq"
  if [[ "$EXTRA_TRIM" -gt 0 ]]; then
    run_cutadapt -u "$EXTRA_TRIM" -u "-$EXTRA_TRIM" --minimum-length "$MIN_LEN_POST" --cores "$threads" -o "$final_fq" "$t5" >/dev/null
  else
    ln -sf "$(basename "$t5")" "$final_fq"
  fi
  final_reads=$(count_fastq_reads "$final_fq")
  final_drop=$((t5_reads - final_reads))
  final_pct=$(awk -v d="$final_drop" -v i="$t5_reads" 'BEGIN { if(i>0) printf "%.2f", (d/i)*100; else print "0.00" }')
  printf "%s\t4_ExtraTrim\t%s\t%s\t%s\t%s%%\n" "$base" "$t5_reads" "$final_reads" "$final_drop" "$final_pct" >> "$run_summary_file"

  # --- STEP 2.5: Post-Trim Length Filtering ---
  final_len_fq="$sample_out/${base}_trim5_u${EXTRA_TRIM}x2_lenFiltered.fastq"
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
  fi


  # --- STEP 3: Initial Alignment ---
  bam_initial="$sample_out/${base}_initial_sorted.bam"
  if [[ "$TRIM_MODE" == "exact" && "${pe_pairs:-0}" -gt 0 ]]; then
    bam_se="$sample_out/${base}_se_sorted.bam"
    bam_pe="$sample_out/${base}_pe_sorted.bam"
    run_log bash -c "{ bwa-mem2 mem $BWA_EXTRA -t $threads '$ref' '$final_fq' | samtools view -Sb - | samtools sort -@ $threads -o '$bam_se'; }"
    run_log bash -c "{ bwa-mem2 mem $BWA_EXTRA -t $threads '$ref' '$pe_r1' '$pe_r2' | samtools view -Sb - | samtools sort -@ $threads -o '$bam_pe'; }"
    run_log bash -c "samtools merge -f -@ $threads '$bam_initial' '$bam_se' '$bam_pe'"
  else
    run_log bash -c "{ bwa-mem2 mem $BWA_EXTRA -t $threads '$ref' '$final_fq' | samtools view -Sb - | samtools sort -@ $threads -o '$bam_initial'; }"
  fi
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
  run_log bash -c "{ $FREEBAYES -f '$ref' -p 1 -E -1 --min-mapping-quality $MAPQ_MIN --min-base-quality $BASEQ_MIN $FB_LIMIT '$bam_baseline' | bcftools norm -m -any --atomize -f '$ref' -Ou 2>/dev/null | bcftools norm -m -any -Ou 2>/dev/null | bcftools view -i 'GT=\"alt\" && ALT!=\"*\"' -Ou | bcftools sort -Ov -o '$vcf_raw_base' 2>/dev/null; }"
  run_log bcftools filter -i "QUAL>$SAFETY_QUAL && INFO/DP>=$MIN_DEPTH" -Ov -o "$vcf_qual_base" "$vcf_raw_base"

  # Baseline Annotation
  if [[ -s "$vcf_qual_base" ]]; then
    sample_bed="$sample_out/${base}_targets.bed"
    # v1.4.1: region BEDs written before v1.4.1 use an older label for the excluded sites; the working copy uses
    # the current one, so only Excluded_Site ever reaches the VCFs and reports
    sed -e "s/PLACEHOLDER_CHROM/$ref_chrom_name/g" -e "s/Blacklist_Site/Excluded_Site/g" "$regions_bed" > "$sample_bed"
    echo '##INFO=<ID=RegionType,Number=1,Type=String,Description="Type of region (HP_Region or Excluded_Site)">' > "$sample_out/hdr.txt"

    vcf_annotated_base="$sample_out/${base}_baseline_annotated_all.vcf"
    run_log bcftools annotate -a "$sample_bed" -c CHROM,FROM,TO,RegionType -h "$sample_out/hdr.txt" -Ov -o "$vcf_annotated_base" "$vcf_qual_base"

    run_log bash -c "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_base' | bcftools filter -e 'POS=7898 || POS=7899 || POS=8595' -Ov -o '$sample_out/${base}_baseline_snps.vcf'"
    run_log bash -c "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_base' | bcftools filter -e 'RegionType=\"HP_Region\" || RegionType=\"Excluded_Site\"' -Ov -o '$sample_out/${base}_baseline_clean.vcf'"
    run_log bcftools filter -i 'RegionType="HP_Region"' -Ov -o "$sample_out/${base}_baseline_homopolymers.vcf" "$vcf_annotated_base"
  fi


  # =========================================================================
  # TRACK B: EXPERIMENTAL (PYTHON TRIMMING/SOFT-CLIPPING)
  # =========================================================================
  echo "[Pipeline] Generating TRIMMED track..." | tee -a "$log_file"
  bam_trimmed="$sample_out/${base}_trimmed_sorted.bam"

  export MIN_ALIGNED_POST MARK_ALIGNED_STATS="$sample_out/${base}_aligned_filter.tsv" MARK_ALIGNED_AMPS="$sample_out/${base}_aligned_filter_by_amplicon.tsv" MARK_PRIMERS_BED="$PRIMERS_BED" MARK_INSERTS_BED="$INSERTS_BED"
  samtools view -h "$bam_initial" | python3 -c '
import os, sys, re

# v1.2.9: the coordinates are read from the BED files, the single source for every primer-dependent step.
# (They were a hard-coded table whose Amp6 and Amp9 product starts were 9 and 6 bases late.)
#   amplicons, tiled_bounds: the callable inserts (CRM_Nested_inserts.bed), 1-based inclusive. Reads are
#     assigned on these spans: assigning on the full product span lets a neighbour with a long primer arm
#     capture reads it then clips away. Adjacent inserts overlap as published for the kit (Vinueza-Espinosa
#     et al. 2023): a base covered by two amplicons is read by two molecules, and both count.
#   products: name, product_start, product_end, forward_primer_end, reverse_primer_start (1-based), from the
#     primer footprints (CRM_Nested_primers_empirical.bed), the fixed read termini seen with this kit.
def read_bed(path):
    rows = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\r\n").split("\t")
        rows[f[3]] = (int(f[1]), int(f[2]))            # 0-based, half-open
    return rows
PB = read_bed(os.environ["MARK_PRIMERS_BED"])
IB = read_bed(os.environ["MARK_INSERTS_BED"])
NAMES = sorted({k.rsplit("_", 1)[0] for k in PB}, key=lambda n: PB[n + "_FWD"][0])
amplicons = [(n, IB[n + "_insert"][0] + 1, IB[n + "_insert"][1]) for n in NAMES]
tiled_bounds = {n: (s, e) for n, s, e in amplicons}
products = [(n, PB[n + "_FWD"][0] + 1, PB[n + "_REV"][1], PB[n + "_FWD"][1], PB[n + "_REV"][0] + 1) for n in NAMES]
PRIMER_TOL = 10
import collections
MIN_ALN = int(os.environ.get("MIN_ALIGNED_POST", "30"))
STATS = os.environ.get("MARK_ALIGNED_STATS", "")
AMPS = os.environ.get("MARK_ALIGNED_AMPS", "")
OUT, bad, famp = [], {}, {}      # famp: fragment -> amplicon, for the per-amplicon counts
def mark(name, flag, reason):
    # v1.2.9: only a primary alignment can remove a read. A supplementary piece (part of a split
    # alignment, e.g. a read with a large deletion) that fails is dropped on its own; the read stays.
    if not flag & 2048:
        bad.setdefault(name, reason)
amp_tag = re.compile(r"_(Amp\d+)")


cigar_re = re.compile(r"(\d+)([MIDNSHPX=])")
# Edge protect is forced to 0 at junctions by using tiled_bounds directly

for line in sys.stdin:
    if line.startswith("@"):
        sys.stdout.write(line)
        continue

    line = line.rstrip("\n")
    parts = line.split("\t")
    if len(parts) < 11: 
        OUT.append((line.split("\t", 1)[0], line + "\n"))
        continue

    flag = int(parts[1])
    # v1.2.8: one read of a pair aligned and the other did not at all - drop both
    if flag & 1 and bool(flag & 4) != bool(flag & 8):
        bad.setdefault(parts[0], "its mate did not align")
    # Skip unmapped and secondary alignments to prevent leaks
    if flag & 4 or flag & 256:
        OUT.append((line.split("\t", 1)[0], line + "\n"))
        continue

    pos = int(parts[3]) - 1 # 0-based
    cigar = parts[5]
    seq = parts[9]
    qual = parts[10]
    
    if cigar == "*" or qual == "*" or seq == "*":
        OUT.append((line.split("\t", 1)[0], line + "\n"))
        continue

    ref_span = sum(int(n) for n, op in cigar_re.findall(cigar) if op in "MDN=X")
    ref_end = pos + ref_span
    if not flag & 2048 and parts[0] not in famp:
        m = amp_tag.search(parts[0])
        famp[parts[0]] = m.group(1) if m else max(amplicons, key=lambda a: min(ref_end, a[2]) - max(pos, a[1]))[0]
    
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
            mark(parts[0], flag, "no template base left inside an amplicon insert")
            continue
        t_start, t_end = tiled_bounds[best_amp]
        keep_start = t_start - 1
        keep_end = t_end

    # off-target guard: the read must still touch at least one amplicon insert
    if not any(min(keep_end, a_end) - max(keep_start, a_start) > 0
               for _, a_start, a_end in amplicons):
        mark(parts[0], flag, "no template base left inside an amplicon insert")
        continue

    if ref_end <= keep_start or pos >= keep_end:
        mark(parts[0], flag, "no template base left inside an amplicon insert")
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
        mark(parts[0], flag, "no template base left inside an amplicon insert")
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

    aligned = sum(int(n) for n, op in cigar_re.findall(final_cigar) if op in "M=X")
    if aligned < MIN_ALN:
        if flag & 2048:
            continue                      # a short supplementary piece: drop the piece, keep the read
        bad.setdefault(parts[0], "fewer than %d template bases aligned" % MIN_ALN)
    OUT.append((parts[0], "\t".join(parts) + "\n"))

# v1.2.7: write everything except reads (and their mates) that failed, then record why
dropped = collections.Counter()
try:
    for q, l in OUT:
        if q in bad:
            dropped[bad[q]] += 1
            continue
        sys.stdout.write(l)
except BrokenPipeError:
    sys.exit(0)
if STATS:
    frags = collections.Counter(bad.values())
    with open(STATS, "w") as fh:
        fh.write("reason\tfragments\talignments\n")
        for k in sorted(set(frags) | set(dropped)):
            fh.write("%s\t%d\t%d\n" % (k, frags.get(k, 0), dropped.get(k, 0)))
if AMPS:
    tot = collections.Counter(famp.values())
    rem = collections.Counter(famp[q] for q in bad if q in famp)
    with open(AMPS, "w") as fh:
        fh.write("amplicon\tfragments\tremoved\tremoved_pct\n")
        for name, _, _ in amplicons:
            if tot[name]:
                fh.write("%s\t%d\t%d\t%.1f\n" % (name, tot[name], rem[name], rem[name] * 100.0 / tot[name]))
' | samtools view -Sb - | samtools sort -n -@ "$threads" - | samtools fixmate -@ "$threads" - - | samtools sort -@ "$threads" -o "$bam_trimmed"

  run_log samtools index "$bam_trimmed"
  aln_drop=$(awk -F'\t' 'NR>1{s+=$2} END{print s+0}' "$MARK_ALIGNED_STATS" 2>/dev/null || echo 0)
  printf "%s\t6_MinAlignedPost(fragments)\t-\t-\t%s\t-\n" "$base" "${aln_drop:-0}" >> "$run_summary_file"
  echo "[Filter] $aln_drop fragments dropped after primer clipping (fewer than $MIN_ALIGNED_POST template bases aligned, none left in an insert, or mate not aligned); see $(basename "$MARK_ALIGNED_STATS")" | tee -a "$log_file"
  flag_aligned_losses "$base" "$MARK_ALIGNED_AMPS"
  flag_low_depth "$base" "$bam_trimmed"
  flag_pileup_cap "$base" "$bam_trimmed"

  # Trimmed Variant Calling
  vcf_raw_trimmed="$sample_out/${base}_trimmed_raw.vcf"
  vcf_qual_trimmed="$sample_out/${base}_trimmed_qual_filtered.vcf"
  run_log bash -c "{ $FREEBAYES -f '$ref' -p 1 -E -1 --min-mapping-quality $MAPQ_MIN --min-base-quality $BASEQ_MIN $FB_LIMIT '$bam_trimmed' | bcftools norm -m -any --atomize -f '$ref' -Ou 2>/dev/null | bcftools norm -m -any -Ou 2>/dev/null | bcftools view -i 'GT=\"alt\" && ALT!=\"*\"' -Ou | bcftools sort -Ov -o '$vcf_raw_trimmed' 2>/dev/null; }"
  run_log bcftools filter -i "QUAL>$SAFETY_QUAL && INFO/DP>=$MIN_DEPTH" -Ov -o "$vcf_qual_trimmed" "$vcf_raw_trimmed"

  # Trimmed Annotation
  if [[ -s "$vcf_qual_trimmed" ]]; then
    vcf_annotated_trimmed="$sample_out/${base}_trimmed_annotated_all.vcf"
    run_log bcftools annotate -a "$sample_bed" -c CHROM,FROM,TO,RegionType -h "$sample_out/hdr.txt" -Ov -o "$vcf_annotated_trimmed" "$vcf_qual_trimmed"

    run_log bash -c "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_trimmed' | bcftools filter -e 'POS=7898 || POS=7899 || POS=8595' -Ov -o '$sample_out/${base}_trimmed_snps.vcf'"
    run_log bash -c "bcftools view -i 'TYPE=\"snp\" && QUAL>$STRICT_QUAL' '$vcf_annotated_trimmed' | bcftools filter -e 'RegionType=\"HP_Region\" || RegionType=\"Excluded_Site\"' -Ov -o '$sample_out/${base}_trimmed_clean.vcf'"
    run_log bcftools filter -i 'RegionType="HP_Region"' -Ov -o "$sample_out/${base}_trimmed_homopolymers.vcf" "$vcf_annotated_trimmed"
    
# =========================================================================
# AUTOMATED HETEROPLASMY EXTRACTION (VAF)
# =========================================================================
echo "[Processing] Extracting heteroplasmy summary for $base..." | tee -a "$log_file"

# Define the output path strictly inside the sample's timestamped folder
het_output="$sample_out/${base}_final_sample_report.tsv"

python3 -c '
import sys, os

# v1.2.8: the sample report is built from the VCFs of the final (primer-clipped) BAM only. The review
# ("baseline") VCFs and the unfiltered *_raw.vcf are not read: they carry primer bases and sites below
# the depth and quality floors, which are not calls.
#   calls                 <base>_trimmed_clean.vcf (SNVs, QUAL > STRICT_QUAL, outside homopolymer regions and excluded sites)
#   review sections       every other record of <base>_trimmed_annotated_all.vcf (QUAL > SAFETY_QUAL, depth >= MIN_DEPTH)
sample_dir, base_name, version, strict_q, safety_q, min_depth = sys.argv[1:7]
het_output = os.path.join(sample_dir, base_name + "_final_sample_report.tsv")


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


def parse_vcf(vcf_path):
    variants = {}
    if not os.path.exists(vcf_path):
        return variants
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            fmt = parts[8].split(":")
            sample_data = parts[9].split(":")
            try:
                tot_dp = int(sample_data[fmt.index("DP")]) if sample_data[fmt.index("DP")] != "." else 0
                ad_parts = sample_data[fmt.index("AD")].split(",")
            except (ValueError, IndexError):
                continue
            if tot_dp == 0:
                continue
            alt_dp = int(ad_parts[1]) if len(ad_parts) > 1 and ad_parts[1] != "." else 0
            vaf = alt_dp * 100.0 / tot_dp
            if not 1.0 <= vaf <= 100.0:
                continue
            info = dict(kv.split("=", 1) for kv in parts[7].split(";") if "=" in kv)
            key = (int(parts[1]), parts[3], parts[4])
            variants[key] = {"lin_pos": int(parts[1]), "ref": parts[3], "alt": parts[4], "vaf": vaf, "dp": tot_dp,
                             "qual": float(parts[5]) if parts[5] != "." else 0.0,
                             "region_type": info.get("RegionType", "")}
    return variants


calls = parse_vcf(os.path.join(sample_dir, base_name + "_trimmed_clean.vcf"))
everything = parse_vcf(os.path.join(sample_dir, base_name + "_trimmed_annotated_all.vcf"))
length_hp, below = {}, {}
for k, v in everything.items():
    if k in calls:
        continue
    snv = len(v["ref"]) == 1 and len(v["alt"]) == 1
    if not snv:
        v["why"] = "length variant; indels are not called"
        length_hp[k] = v
    elif v["region_type"] == "HP_Region":
        v["why"] = "homopolymer region; not called"
        length_hp[k] = v
    elif v["region_type"] == "Excluded_Site":
        v["why"] = "excluded site; not called"
        below[k] = v
    else:
        v["why"] = "QUAL %.0f, calls need > %s" % (v["qual"], strict_q)
        below[k] = v

HEADER = ("Region / Locus\tRange Covered\trCRS Reference Nucleotide\tVariant Called\tHeteroplasmy % (VAF)\t"
          "Notes / HGVS Equivalent\tRead Depth\n")

with open(het_output, "w") as out:
    def write_variants(var_dict):
        for v in sorted(var_dict.values(), key=lambda x: (x["lin_pos"] + 8284) % 16569 or 16569):
            lin_pos, ref, alt, vaf = v["lin_pos"], v["ref"], v["alt"], v["vaf"]
            real_pos = (lin_pos + 8284) % 16569
            if real_pos == 0: real_pos = 16569
            region_name, range_cov = get_region(real_pos)
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
                    variant_called = f"m.{del_start}_{del_start + len(del_seq) - 1}del"
                notes = "Deletion"
                het_status = f"Homoplasmic {vaf:.0f}%" if vaf > 95.0 else f"{vaf:.0f}% LHP"
            else:
                disp_ref = "-"
                variant_called = f"m.{real_pos}.1{alt[len(ref):]}"
                notes = "Insertion"
                het_status = f"Homoplasmic {vaf:.0f}%" if vaf > 95.0 else f"{vaf:.0f}% LHP"
            if "why" in v:
                notes = notes + "; " + v["why"]
            depth = v["dp"]
            out.write(f"{region_name}\t{range_cov}\t{disp_ref}\t{variant_called}\t{het_status}\t{notes}\t{depth}\n")

    out.write(f"# MARK v{version} sample report: {base_name}\n")
    out.write(f"# From the final (primer-clipped) BAM only. Calls = {base_name}_trimmed_clean.vcf: SNVs with QUAL > {strict_q} "
              f"and depth >= {min_depth}, outside homopolymer regions and excluded sites. The review sections list the other "
              f"sites with QUAL > {safety_q} ({base_name}_trimmed_annotated_all.vcf); they are not calls. "
              f"Heteroplasmy % = reads with the variant / read depth.\n")
    out.write("### MAIN RESULTS (the calls) ###\n")
    out.write(HEADER)
    write_variants(calls)
    out.write("\n### LENGTH VARIANTS AND HOMOPOLYMER SITES (for review, not calls) ###\n")
    out.write(HEADER)
    write_variants(length_hp)
    out.write("\n### BELOW THE CALL THRESHOLD OR AT AN EXCLUDED SITE (for review, not calls) ###\n")
    out.write(HEADER)
    write_variants(below)
' "$sample_out" "$base" "$VERSION" "$STRICT_QUAL" "$SAFETY_QUAL" "$MIN_DEPTH"

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

  run_log fastqc "$bam_trimmed" -o "$sample_out" > /dev/null 2>&1
  echo "Completed Diagnostic Split: $base" | tee -a "$log_file"
done

if ! $found_any; then
  echo "No paired FASTQs found for '$input_path'."
  exit 0
fi

if [[ -f "$run_out/WARNINGS.txt" ]]; then
  {
    echo "################################################################################"
    echo "## THIS RUN HAS WARNINGS - see $run_out/WARNINGS.txt"
    awk -F'\t' 'NR>2 { printf "##   %s %s: %s (%s, flag level %s)\n", $1, ($3 == "all" ? "" : $3), $2, $4, $5 }' "$run_out/WARNINGS.txt"
    echo "################################################################################"
  } >&2
fi
echo "All outputs written to: $run_out"
