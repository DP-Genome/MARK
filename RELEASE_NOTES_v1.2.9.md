# MARK v1.2.9 — release notes

Everything here is new since **v1.2.6**, the previous release. It covers v1.2.7, v1.2.8 and v1.2.9,
which were developed and validated together and not released separately.

**The variant-calling settings and filters are unchanged.** v1.2.9 changes three things:
- which Illumina reads are kept;
- which aligned reads reach the final BAM;
- what a run tells you when something needs attention.

The variant calls were identical to v1.2.6 in all 61 validation libraries.

---

## 1. Illumina keeps only read pairs that merge

A fragment is kept only when read 1 and read 2 both pass filtering **and** merge into one fragment.

| | v1.2.6 | v1.2.9 |
|---|---|---|
| exact mode (default) | merged reads, plus unmerged pairs whose two reads start with the primer pair of one amplicon | **merged reads only**; pairs that merge in neither merge pass are counted and dropped |
| standard mode | merged reads, plus single reads whose mate failed filtering (mate rescue) | **merged reads only**; when one read fails, both go |

**Why:**
- Every read in the Illumina final BAM is now one whole amplicon read from both ends.
- The unmerged pairs were mostly Amp10 pairs whose read 2 loses phase after the rCRS 431–438 C-tract.
  The degraded part was soft-clipped out of the calls, but the bases just before it carried extra
  errors that did reach the caller.

**What to expect:**
- **Depth falls in Amp10, and only there.** Its 218 bp product needs 84 bases of overlap from two
  151-base reads, so many Amp10 pairs never merge. In the validation set Amp10 kept 9–57% of its depth
  and stayed at 300x or more, except one library (below). The other nine amplicons keep 85–100%.
- **A library whose read 2 fails across Amp10 can lose that amplicon.** In the validation set this
  happened to one library, DNA007, whose Amp10 fell to 1–2x. The new depth warning (section 3) names it.
- **The Illumina background is lower.** In the noisiest 1% of sites, non-reference bases fell from
  1.58% to 0.40% in one library and from 1.72% to 0.97% in another.

The run summary records the dropped pairs as `2b_UnmergedPairs_dropped`.

## 2. Reads that align almost nothing are removed after primer clipping

After the primer bases are soft-clipped, **a read must still have at least 30 template bases aligned**
(`MIN_ALIGNED_POST`, default 30), or it is removed. This catches reads that pass every length filter
but hold almost no real sequence, for example an Illumina read that turns to poly-G right after its
primer (`27S3M121S`).
- When the clipping step leaves a read with no template base inside an amplicon insert, it is removed,
  as before.
- Each library writes `<library>_aligned_filter.tsv` (removals by reason) and
  `<library>_aligned_filter_by_amplicon.tsv` (per amplicon).
- The run summary gains a `6_MinAlignedPost` row.
- The review ("baseline") BAM is not filtered.
- **The clipping step reads the amplicon coordinates from the BED files** (`PRIMERS_BED` and
  `INSERTS_BED`), the same files the rest of the pipeline uses, instead of a table written into the
  scripts. Both files must be present. A custom primer or insert BED now applies to the clipping as well.
- **Only a read's primary alignment can remove it.** A supplementary alignment (one piece of a split
  read) is never checked on its own, so a read whose large deletion makes the aligner split it is kept.

Real reads align 40 or more template bases. In the validation set this check removed 0–26 nanopore
reads per run and model, and a few Illumina reads per library.

## 3. Warnings: an amplicon below the calling depth, or losing reads after alignment

`WARNINGS.txt` at the top of the run folder now has three checks. None of them stops the run.

| check | flagged when | setting |
|---|---|---|
| primer-pair share (unchanged) | fewer than 25% (Illumina) or 10% (nanopore) of reads carry a primer pair of one amplicon | `PRIMER_WARN_PCT` |
| **amplicon depth (new)** | any insert position of an amplicon is below `MIN_DEPTH` (10x) in the final BAM, so variants there cannot be called | `MIN_DEPTH` |
| **reads removed after alignment (new)** | the check in section 2 removes 20% or more of an amplicon's fragments, and at least 10 | `ALIGNED_WARN_PCT`, `ALIGNED_WARN_MIN` |

- **Deletions do not trigger the depth check.** A read with a deletion at a position counts as covering
  it, so a homozygous deletion is not mistaken for a gap.
- **Unclassified barcode bins are skipped**, because they are not samples.
- **Every library writes `<library>_amplicon_depth.tsv`**, with the mean and minimum depth per amplicon.

**`WARNINGS.txt` has a new column layout**, the same for every check:

```
library  check  amplicon  value  flag_level  count  of  detail
```

The values carry their unit, for example `96.4%` or `1x`. A script that read the v1.2.6 columns
(`primer_pair_pct`, `warn_below_pct`, ...) needs updating.

## 4. The per-library sample report reports the final calls

`<library>_final_sample_report.tsv` used to merge every VCF in the library folder, including the
review ("baseline") and unfiltered (`*_raw.vcf`) call sets. As a result it could list artefacts that
are not calls. One example is `m.15985`, one base before the Amp1 primer, which came from the unclipped
review alignment.

It is now built from the final alignment only, in three sections:
1. **the calls**: exactly the `<library>_trimmed_clean.vcf` calls, with heteroplasmy and read depth;
2. **length variants and homopolymer sites**, for review, each with the reason it is not a call;
3. **sites below the call threshold or at a blacklisted site**, for review.

A read-depth column is added. The run summary header now records `MIN_ALIGNED_POST`, `ALIGNED_WARN`
and `DEPTH_WARN`.

## 5. Length settings: the same on both platforms where possible

| | v1.2.6 | v1.2.9 |
|---|---|---|
| Illumina `MAX_LEN` | a fastp pass removing pairs longer than 1,500 bp | **removed**: reads are 2 × 151 bp, so a maximum cannot apply. Empty reads, which this pass used to remove, now fail the 90-base floor at merging |
| Illumina, exact mode | no product-length window | **products of 90–300 bp**, the same window as nanopore (run summary row `2d_ProductLength`) |
| Illumina `MIN_LEN_POST` | 30 | **90**, as on nanopore |
| nanopore | `MAX_LEN` 1,500; products 90–300 bp | unchanged |

- **Exact mode:** every product kept is one whole amplicon, 144–237 bp, so the 90–300 bp window removes
  nothing.
- **Standard mode:** the new 90 bp floor removed 0.25–0.85% of fragments in the validation set, with no
  change to the calls.

## 6. Launcher (`MARKLaunch.py`)

- **`MIN_ALIGNED_POST`** is in the parameter table, with a "?" help card.
- **Help cards:**
  - `MAX_LEN` is marked nanopore-only, and `MIN_LEN_POST` and `MAX_LEN_POST` apply to both platforms.
  - The `TRIM_MODE` card says that exact mode is the default.
  - The `MIN_DEPTH` card mentions the depth warning.
- **The final pop-up reports warnings.** When the run folder has a `WARNINGS.txt`, "Run pipeline" and
  the auto cycle end with a warning pop-up titled "finished with N warning(s)". It lists each warning
  and gives the paths to `WARNINGS.txt` and the results.
- **`Final_Pipeline_Results/` says what the result is.** `corrected_vcfs/` and `cleaned_bams/` keep
  their names and file names, so `VCF_OrganizerApp` works as before.

```
Final_Pipeline_Results/
  README.txt              what every folder holds
  WARNINGS.txt            when the run has warnings
  run_summary.txt
  corrected_vcfs/         the calls only: <library>_trimmed_clean_corrected.vcf, _trimmed_snps_corrected.vcf
  cleaned_bams/           the final BAM only: <library>_trimmed_sorted_cleaned.bam (+ .bai)
  sample_reports/         <library>_final_sample_report.tsv
  qc/                     <library>_amplicon_depth.tsv, <library>_aligned_filter_by_amplicon.tsv
  reference/              linearized_mtdna.fasta (+ .fai), to open the BAMs in IGV
  review/corrected_vcfs/  every other VCF (review-BAM calls, raw and intermediate call sets, homopolymers)
  review/cleaned_bams/    the review ("baseline") BAM
```

The alignment intermediates (`initial`, `pe` and `se` BAMs) stay in `sorted_bams/`.

- **The corrected VCFs are fully in rCRS coordinates:**
  - `CHROM` and the `##contig` line are the original reference's name (`NC_012920.1`, 16,569 bp);
  - `##reference` points to the original reference;
  - the records are sorted by position.
  
  They pass `bcftools norm --check-ref` against `rCRS.fasta`, and can be bgzipped and indexed. Corrected
  VCFs made by earlier versions say `NC_012920.1_linearized`, so compare old and new files by position.

## 7. New and changed environment variables

| variable | default | pipeline | what it does |
|---|---|---|---|
| `MIN_ALIGNED_POST` | `30` | both | template bases a read must align after primer clipping (section 2) |
| `ALIGNED_WARN_PCT` | `20` | both | flag an amplicon losing this % of its fragments to that check |
| `ALIGNED_WARN_MIN` | `10` | both | ... and at least this many fragments |
| `INSERTS_BED` | `CRM_Nested_inserts.bed` | both | amplicon inserts, for the clipping step and the depth check |
| `MIN_LEN_POST` | `90` (Illumina was `30`) | both | lower end of the product-length window |
| `MAX_LEN` | `1500` | **nanopore only** (removed from Illumina) | raw read length limit |

All other variables keep their v1.2.6 names and defaults.

## 8. Validation

v1.2.9 was run on 61 libraries (4 Illumina, 57 nanopore, with the Flongle run left out) through the
launcher's auto cycle with every parameter at its default, and compared with v1.2.6, v1.2.7 and v1.2.8:
- **Calls:** identical in all 61 libraries, in both the SNP and the clean VCF.
- **Nanopore:** final BAMs identical to v1.2.8, read for read.
- **Illumina:** Amp10 depth as in section 1.
- **BED coordinates and the primary-alignment rule (section 2)** were checked again on all 61 libraries:
  calls identical in all 61; 17 more Illumina fragments kept, all split alignments (13 of them NIST B reads
  carrying a deletion), none lost; 33 nanopore reads clipped a few bases differently.
- **Warnings:** the run produced exactly one, `DNA007 Illumina Amp10: amplicon depth below MIN_DEPTH`
  (156 of 164 insert positions at 1–2x).
- **Tool versions:** fastp 1.3.7, fastplong 0.7.1, cutadapt 5.2, bwa-mem2 2.3, minimap2 2.31,
  samtools 1.24, bcftools 1.24, FastQC 0.12.1, as the conda recipe requires. The bioconda bwa-mem2 2.3
  binary reports its own version as "2.2.1", so BAM headers (`@PG ... VN:2.2.1`) and `bwa-mem2 version`
  show 2.2.1 even though the 2.3 package is installed.

## 9. Known limitations

- **Illumina depth over Amp10 depends on read 2 surviving the rCRS 431–438 C-tract.** A library whose
  read 2 fails there loses Amp10. The depth warning names it, and the nanopore data covers the region.
- **No mate rescue on either Illumina path.** Both reads must pass and merge.
- **Reported calls are SNVs.** Indels and length variants are listed in the sample report for review,
  not called, as in earlier releases.

## 10. Upgrading

- **Calling:** nothing to change; the variant filters are the same.
- **Illumina:** expect lower Amp10 depth, and check `WARNINGS.txt` for any amplicon below the calling
  depth.
- **Scripts that read MARK outputs:** update any that read `WARNINGS.txt` (new columns),
  `final_sample_report.tsv` (new sections and depth column), the corrected VCFs (new contig name), or
  that expect intermediate BAMs or every VCF directly in `Final_Pipeline_Results/`.
