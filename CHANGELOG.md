# MARK changelog

What changed in each version of MARK, newest first. The notes for each version are also on the
[GitHub Releases page](https://github.com/DP-Genome/MARK/releases).

- [v1.7](#v17) (2026-11-09)
- [v1.6.2](#v162) (2026-10-29)
- [v1.6.1](#v161) (2026-10-27)
- [v1.6](#v16) (2026-10-27)
- [v1.5](#v15) (2026-10-22)
- [v1.4.4](#v144) (2026-10-22)
- [v1.4.3](#v143) (2026-10-16)
- [v1.4.2](#v142) (2026-10-16)
- [v1.4.1](#v141) (2026-10-15)
- [v1.4](#v14) (2026-10-14)
- [v1.2.10](#v1210)
- [v1.2.9](#v129)
- [v1.2.6](#v126) (2026-09-23)
- [v1.2.5](#v125)
- [v1.2.4](#v124)
- [v1.2.3](#v123)
- [v1.2.1](#v121) (2026-09-15)

---

<a id="v17"></a>

## MARK v1.7: mixture check, heteroplasmy report fix

v1.7 changes four things relative to v1.6.2. The main variant calls (the consensus profile) are unchanged in all reference datasets.

1. **A mixture check.** MARK now looks for evidence of a mixture: a second mtDNA haplotype in the library, from two people's DNA or from contamination.
2. **A fix to the heteroplasmy report.** freebayes reports one combined allele when a homopolymer length change sits next to a base change. MARK split these alleles base by base (`bcftools norm --atomize`). That shifted every base after the indel, which created substitutions no read carries. It also split one base's reads across records.
   - Examples on nanopore: NIST B/C 71T/72G at about 7%, and 73G at about 90% (homoplasmic in the reads).
   - Example on nanopore: DNA 007 460C at about 70% (homoplasmic C8).
   - Example in a simulated test: a 30% heteroplasmy next to NIST C's AC-repeat deletion was missed.

   The new `MARK_het_decompose.py` splits each allele by its CIGAR. A substitution's fraction now counts only the reads that have a base at that position; reads with a deletion there are a length variant. Only the heteroplasmy report changes.
3. **Completed homopolymer mask** (`linearized_regions.bed`). Five homopolymers of 5 or more bases in the panel were not marked as homopolymer regions: G6 at 66–71, A5 at 357–361, C6 at 433–438, C6 at 16375–16380 and G5 at 16470–16474. Each is now marked, from the base before it, like the others.
   - Effect: nanopore 1-base deletions in these tracts, listed at 5–12%, now carry the homopolymer-region label (for review).
   - No variant call in the 60 reference datasets lies in these tracts.
4. **A clearer warning when only one amplicon is involved.** When all the linked bases are in one amplicon, the warning says "one amplicon only". That pattern can be contamination concentrated in that amplicon, or an artefact of it.
   - In the reference run this applies to the three DNA 007 Illumina libraries: amplicon 10, the Illumina-only background.

### How it works
A second person's bases at the positions where the two haplotypes differ sit **on the same reads**. Sequencing errors and DNA damage also leave minor bases, but each sits on its own read. The new script `MARK_mixture_check.py` runs on the final, primer-clipped BAM:

1. A candidate is the second most common base at a position, outside homopolymer regions and excluded sites, carried by at least 1% of reads and at least 5 reads.
2. Two candidates are **linked** when they meet all of these:
   - at least 50 reads cover both positions;
   - at least 5 of those reads carry both minor bases;
   - that count is at least 70% of the smaller minor-base count;
   - the share of covering reads that carry both is at least `MIX_MIN_AF`.
3. One linked pair is enough to flag the library. The estimated share of the second haplotype is the median over the linked pairs.

### Outputs
- **WARNINGS.txt:** "evidence of a mixture (second haplotype)", with the amplicons, the estimated share and the linked bases.
- **The sample report:** a new section "MIXTURE CHECK: A SECOND HAPLOTYPE?". It is written for every library, flagged or not.
- **`<lib>_mixture_check.tsv`:** the linked pairs and their read counts, copied to `sample_reports/`.
- **Run summary and launcher:** a run-summary line, and the console banner during the run.

### Settings
- `MIX_MIN_AF` (default 0.02): the smallest share of the second haplotype that is reported. `0` turns the check off.
- `HET_MIN_AF` (unchanged, default 0.05) can still be raised or lowered.

### Validation
Reads of 2800M, DNA 007, NIST B and NIST C were mixed at the FASTQ level, with every ordered pair at 0.5, 1, 2, 5 and 10%, and DNA 007/2800M also at 20 and 50%. That gives 63 mixtures per platform, plus each sample alone. Every mixture went through full MARK v1.7.

**Mixtures flagged by the check**

| Second sample | Illumina | Nanopore |
|---|---|---|
| none (sample alone) | 0/4 | 0/4 |
| 0.5% | 3/12 | 1/12 |
| 1% | 4/12 | 6/12 |
| 2% | 10/12 | 8/12 |
| 5% | 12/12 | 12/12 |
| 10% | 12/12 | 12/12 |

**Estimated share in flagged libraries (median):**
- 2%: 2.9% on Illumina, 2.5% on nanopore;
- 5%: 5.1% on both;
- 10%: 10.0% on Illumina, 10.6% on nanopore.

The share can differ between amplicons, because each amplicon's depth differs between the two samples.

**HET_MIN_AF on the DNA 007 + 2800M mixtures (rCRS 460 left out):**
- At 5% (default), every site with a true fraction of 5% or more was reported, and none below 5%.
- At 2% and at 10%, the cut-off moved with the setting. The only exceptions were 8 nanopore sites whose true fraction lay just under the setting but were measured just over it:
  - 5 at 1.7–1.98% true, measured 2.0–2.5%;
  - 3 at 9.6–9.9% true, measured 10.0–10.3%.

**Heteroplasmy report, against the reads (the share of reads with the second base, counted from the final alignment at base and mapping quality ≥ 20):**
- **Real data (60 datasets).** Expected means the reads show 5–95%.
  - Illumina: 29 of 29 expected were reported outside homopolymer regions, with 0 reported beyond them.
  - Nanopore: 22 of 24 were reported, with 0 beyond them. The 2 misses are 16373C, a run artefact, at a borderline 5.3–5.5%.
- **Simulated heteroplasmy.** A second base was written into a known share of real reads: 10 levels from 1% to 50%, one position per amplicon, in each of the four samples, on both platforms.
  - At HET_MIN_AF 2%, 5% and 10%, every one at or above the setting was reported and none below it: 240 of 240 as expected.
- **DNA 007 + 2800M mixtures, at 5%.**
  - Illumina: 54 of 54 positions at ≥5% were reported, and 0 of 63 below 5%.
  - Nanopore: 53 of 53 at ≥5% were reported, and 0 of 64 below 5%.

**Mixture warning, scored against the true share in the reads.** For each mixture, the true share is measured on the reads covering two positions where the samples differ (the reads the check uses), taking the highest such pair.
- At 2% or more: Illumina 44 of 45 warned, nanopore 42 of 42.
- The one miss is in Illumina amplicon 10. There, DNA 007's own reads at 460 are already a C/T mix.
- Below 2%: no warnings on either platform (Illumina 0 of 18, nanopore 0 of 21).
- Libraries alone: 0 of 8.
- MARK's estimated share matches the true share on the same reads, except in four Illumina mixtures with DNA 007 (its amplicon 10).

**Real data (MARK v1.6.2 reference run, 60 datasets):** the check flagged only the three DNA 007 Illumina libraries. These carry an Illumina-only linked background haplotype in amplicon 10 at about 4%, which is under separate investigation.

---

<a id="v162"></a>

## MARK v1.6.2

**Primer bases are always hard-clipped.** The `CLIP_MODE` setting is removed: after alignment, the primer bases at
each product end are removed from the read (CIGAR `H`), which has been the default since v1.4.3. The soft-clip
option, its launcher drop-down and its help text are gone, and the run summary states the clipping as
"PRIMER CLIPPING: primer bases hard-clipped". A `CLIP_MODE` value left in the environment is ignored.

The pipelines, the calls and the heteroplasmy report are otherwise identical to v1.6.1: re-running all 25 reference
runs (63 libraries) gave the same alignments, VCFs, read counts and trimmed reads, and sample reports that differ only
in their title line ("MARK v1.6.2 sample report").

---

<a id="v161"></a>

## MARK v1.6.1

**Launcher layout only.** The Activity Log strip is 132 px tall instead of 156 px, so the whole parameter table,
including the new HET_MIN_AF row and the last row (FREEBAYES), is visible without scrolling when the window is
maximized. The pipelines, the calls and the heteroplasmy report are identical to v1.6.

---

<a id="v16"></a>

## MARK v1.6

**Heteroplasmy report.** Every allele at 5% or more of the reads is listed, with its fraction, in a VCF and in the
sample report. The calls are identical to v1.5.

### What changed
- **Why a second freebayes run:** the calls come from a one-copy model, which reports the majority base at each
  position. A 10% heteroplasmy at a non-variant site, or the minor allele of a 70/30 site, is never shown.
- **What v1.6 runs:** freebayes again on the final (primer-clipped) BAM, in pooled-continuous mode. This mode gives
  allele fractions, not genotypes.
- **What is reported:** every allele that meets all three of these:
  - fraction of at least `HET_MIN_AF`;
  - fraction of at most 1 − `HET_MIN_AF`;
  - depth of `MIN_DEPTH` or more.
- **Where it is reported:**
  - **`<library>_trimmed_heteroplasmy.vcf`:** the fraction is in `INFO/HET_AF`, and homopolymer regions and excluded
    sites are annotated (`RegionType`). After post-processing it is in `corrected_vcfs/` (rCRS).
  - **The sample report:** a new section, `### HETEROPLASMY: ALLELES AT 5% OR MORE (for review, not calls) ###`. Each
    entry is point (PHP) or length (LHP) heteroplasmy, in forensic notation, with homopolymer regions and excluded
    sites marked.
- **Complex alleles:** freebayes sometimes writes an insertion next to a base change as one allele. Its single-base
  parts are added up, so the fraction is the share of all reads at the position that carry the base.
- **New setting `HET_MIN_AF`:** default `0.05`; `0` turns the report off. It is in the launcher's parameter table
  (with help text), `--help`, the README settings table and the run summary.
- **freebayes requirement:** the report needs freebayes. With `VARIANT_CALLER=bcftools` and no freebayes installed,
  the report is skipped with a note in the log.

### Unchanged
Trimming, alignment, calls (`_snps`, `_clean`), filters, and the other outputs are the same as v1.5.

### Validation
Investigation: `analysis_and_comparisons/v1.5/investigations/heteroplasmy_indels_20261027`.
Test runs: `analysis_and_comparisons/v1.6/tests/heteroplasmy_20261027`.

**In-silico mixtures** of 2800M and DNA 007 reads (2–50%, both directions, Illumina and nanopore, 18 mixtures):

| Heteroplasmies | Reported |
|---|---|
| True fraction ≥ 5% (110) | 109 (the one missed was at 5.0%); fractions within 0.8 percentage points |
| True fraction 2% | none |

Mutserve 2.0.3 gives the same results, but only on alignments with reads on both strands. MARK's alignments have
every read on one strand, and on them Mutserve reports no heteroplasmy at all. LoFreq and Mutect2 are less accurate.

**Real libraries:**
- **2800M:** no heteroplasmy on either platform. The one exception is 16373C at 5–7%, which appears in every sample of
  the 2025 MinION run (a run artefact).
- **DNA 007:**
  - 385G at 5–8% in every library on both platforms. This is a probable low-level heteroplasmy; other samples are
    at ≤ 0.3%.
  - 460C at about 45% (Illumina) and about 70% (nanopore), from the C8 tract.
  - The 573 C-tract length heteroplasmy (+4C dominant).
  - Amp10 read-2 noise at 437–498 on Illumina.

**End-to-end runs:** command line, then `--postprocess`.
- **Nanopore:** 2800M and DNA 007, SUP v5.2 (2025), 6 libraries. Exit 0.
- **Illumina:** DNA007-R1 and 2800M. Exit 0.
- **Outputs:** heteroplasmy VCFs are in `corrected_vcfs/`, and the sample reports have the new section. The values
  match the standalone analysis.

---

<a id="v15"></a>

## MARK v1.5

**Choose the variant caller: freebayes (the default), bcftools, or both, with a per-library comparison of the two.
With the default, the calls are identical to v1.4.4.**

### What changed
- **New setting, `VARIANT_CALLER`** (the first row of the launcher, a drop-down; or `VARIANT_CALLER=... MARK.sh`):
  - **`freebayes` (default):** as v1.4 to v1.4.4. freebayes realigns the reads over each candidate region and judges
    whole haplotypes.
  - **`bcftools`:** `bcftools mpileup` + `call`, as MARK used up to v1.2.10, through the same filters and with the same
    file names. freebayes does not need to be installed for this mode.
  - **`both`:** freebayes makes the calls (the usual VCFs and the sample report). bcftools also runs on the same final
    BAM, through the same filters, and MARK compares the two.
- **What `both` writes, per library:**
  - `<library>_trimmed_bcftools_snps.vcf` and `_clean.vcf`: bcftools' own call sets;
  - `<library>_trimmed_snps_both_callers.vcf`: every SNV in either call set, tagged `CALLERS=both`, `freebayes_only`
    or `bcftools_only`;
  - `<library>_caller_comparison.txt`: a plain-text table of every such variant with each caller's verdict
    (called / below threshold / not found), QUAL and depth, and a note to review the alignment where only one caller
    called it.
  - Post-processing converts the new VCFs to rCRS coordinates with the rest (`corrected_vcfs/`) and files the
    comparison with the sample reports (`sample_reports/`).
- **Sample report and run summary** name the caller that made the calls.
- **Command line, end to end:** `MARKLaunch.py --postprocess RUN_FOLDER` runs the four post-processing steps (VCFs to
  rCRS, cleaned BAMs, `Final_Pipeline_Results/`) without the GUI, so a whole analysis can be scripted:
  ```bash
  VARIANT_CALLER=both RUN_NAME=run1 OUTPUT_DIR=/data/out MARK.sh /data/fastqs
  MARKLaunch.py --postprocess /data/out/run1
  ```
- **`--help`:** `MARK.sh --help` and `MARK-I.sh --help` now list every setting with its default, grouped (calling,
  trimming, read filters, warnings, files), and exit with status 0. `MARKLaunch.py --help` describes the launcher's
  options.
- **README:** a command-line guide and a settings reference table covering every parameter, the same settings as the
  launcher's parameter table.
- **Launcher:** the compact Activity Log of v1.4.4, sized so the longer parameter table (one more row) still fits
  without scrolling when the window is maximized.

### Thresholds
Both callers go through the same filters (depth ≥ `MIN_DEPTH`, QUAL > `SAFETY_QUAL` for review, SNVs with QUAL >
`STRICT_QUAL` for calls). These were calibrated for freebayes; bcftools' QUAL is on a different scale (it tops out
around 228), so a bcftools-only or bcftools-missed variant near the threshold should be read in the alignment.

### Validation
- **Default unchanged:** with `VARIANT_CALLER=freebayes`, the calls are identical to the v1.4 reference run.
- **DNA 007 · Illumina 2026 (3 replicates), `both`:** 11 of 14 SNVs are called by both in each; 460C (in the 456–463 C tract) is `freebayes_only` in every
  replicate; 572A and 573T are `bcftools_only` (C-tract artefacts inside a homopolymer region).
- **DNA 007 · MinION 2026 SUP (3 libraries), `both`:** in each, 11 of 12 SNVs are called by both; 16223T is
  `freebayes_only` (bcftools found it, but with QUAL below `STRICT_QUAL`).
- **`bcftools` alone:** runs end to end on the DNA 007 nanopore libraries; the freebayes check is skipped in this mode.
- **Release check (this build), DNA 007 subsets:** command line, pipeline then `--postprocess`, for nanopore with each
  of `freebayes`, `bcftools` and `both`, and Illumina with `both` (one R1 file) and `bcftools` (a folder of pairs): all
  finished with exit status 0 and complete `Final_Pipeline_Results/`. Each caller's calls are identical whether it runs
  alone or in `both`. The launcher's AUTO cycle with `VARIANT_CALLER=both` set in its parameter table produced the
  same set of outputs.

---

<a id="v144"></a>

## MARK v1.4.4

**Launcher layout only: the Activity Log is a compact strip, so the whole parameter table fits on screen. The pipeline
and its results are identical to v1.4.3.**

### What changed
- **Activity Log:** a fixed strip about 170 px tall (it scrolls) at the foot of the dashboard, instead of a panel that
  stretched with the window and took about a quarter of it (290 px on a maximized MacBook Pro window).
- **Parameter table:** the tabs take all the spare height. With the window maximized, every parameter row, down to
  FREEBAYES, is visible without scrolling.
- **README:** the dashboard screenshot shows the current launcher (it was from v1.1.5).
- **Version:** 1.4.4 in MARK.sh, MARK-I.sh, the launcher title and header, and the conda recipe.

### Unchanged
MARK.sh and MARK-I.sh differ from v1.4.3 only in the version string. Trimming, clipping, variant calling, filters,
reports and output files are the same.

### Validation
- **Layout:** on a maximized 2056 × 1188 pt window, the Activity Log is 168 pt and the parameter table needs 823 of the
  836 pt available, so it does not scroll (v1.4.3: log 290 pt, table area 765 pt, table scrolled).
- **Pipeline:** `bash -n` passes for MARK.sh and MARK-I.sh; their diff against v1.4.3 is the version line.

---

<a id="v143"></a>

## MARK v1.4.3

**Primer bases are hard-clipped after alignment by default. The variant calls are identical to v1.4.2.**

### What changed
- **New setting, `CLIP_MODE`:** controls how MARK clips the primer bases at each product end after alignment.
  - **`hard` (the new default):** the primer bases are removed from the read (CIGAR `H`; SEQ and QUAL shortened).
  - **`soft`:** the primer bases stay in the BAM as soft clips (CIGAR `S`), as up to v1.4.2.
  - Either way only insert bases reach variant calling, and the calls are the same.
- **Hard clips from the aligner:** any hard clips the aligner made (supplementary pieces of split alignments) are
  kept and added to the clip length.
- **Launcher:** `CLIP_MODE` is a drop-down (hard / soft) with help text.
- **Run summary:** records the clip mode.
- **Documentation:** the README and the launcher's Final_Pipeline_Results README describe the final BAMs as
  hard-clipped.

### Why hard clipping
- The final (`trimmed`/`cleaned`) BAMs are about 20% smaller.
- They hold only the bases that variant calling uses.
- The untrimmed baseline BAM still shows the primers for review.

### Validation
- **Hard vs soft:** on v1.4.1 with `CLIP_MODE=hard`, compared with soft clipping across 10 libraries (DNA 007 ·
  Illumina 2026 ×3, DNA 007 · MinION 2025 SUP v5.2 ×3, NIST C · MinION 2023 SUP v5.2 ×3, NIST B · Illumina 2025):
  - all 240 VCFs had identical records;
  - the trimmed BAMs had the same reads and the same aligned bases;
  - every soft clip became a hard clip.
- **v1.4.3 with default settings:** DNA 007 · MinION 2025 SUP v5.2 (3 libraries) and 2800M · Illumina 2025. All
  192 VCFs are identical to v1.4.2, the final BAMs carry hard clips only, and the run summary records CLIP_MODE=hard.
- **`CLIP_MODE=soft`:** reproduces the v1.4.2 CIGARs (e.g. 22S110M23S for a whole Amp2 product; hard: 22H110M23H).

---

<a id="v142"></a>

## MARK v1.4.2

**Sample reports in forensic notation. The variant calls and VCFs are identical to v1.4.1.**

### What changed
The sample report (`Final_Pipeline_Results/sample_reports/<library>_final_sample_report.tsv`) writes variants the way
forensic mtDNA profiles are written (ISFG and EMPOP conventions):
- **SNVs:** `152C` (was `m.152C`). The HGVS form, `m.152T>C`, stays in the Notes column.
- **Insertions:** moved to the 3-prime end of the repeat they lie in and numbered there. The VCF record 310 T>TC
  becomes `315.1C`, and a three-C insertion in the 568–573 tract becomes `573.1C 573.2C 573.3C`.
- **Deletions:** also moved to the 3-prime end, one entry per deleted position. A CA deletion in the 514–523 repeat
  becomes `523del 524del`.
- **Traceability:** the Notes column keeps the VCF record of each insertion and deletion, e.g.
  `Insertion (VCF: 310 T>TC)`.

The VCFs are unchanged: they keep the standard left-aligned representation. Insertions and deletions remain
review-only (the "LENGTH VARIANTS AND HOMOPOLYMER SITES" section); the call set is still SNVs only.

### Known limitation
Around the HV2 C-stretch (303–315), freebayes can split one complex haplotype into two records, e.g. 310 T>C plus
an insertion. Each record is written on its own, which can give an unusual notation such as `310.1T 310.2C …`.
These sites are review-only; read them against the alignment, or align the haplotype phylogenetically (e.g. EMPOP).

### Validation
- **Unit test:** VCF records to notation for 310 T>TC → 315.1C, 152 T>C → 152C, +CCC at 567 → 573.1C–573.3C, GCA>G
  at 513 → 523del 524del, and C insertions in the HV1 C-tract → 16188.1C and 16193.1C.
- **Reference run:** reports rebuilt from the reference-run VCFs for 2800M (152C, 263G, 477C, 16519C; 315.1C under
  review), DNA 007 nanopore (573.1C–573.4C under review) and NIST C (523del 524del under review).
- **Full pipeline test** against the previous version, on DNA 007 · MinION 2025 SUP v5.2 (3 libraries), 2800M ·
  Illumina 2025 and NIST C · Illumina 2025: all 240 VCFs have identical records.

---

<a id="v141"></a>

## MARK v1.4.1

**Terminology change only. The variant calls are identical to v1.4.**

### What changed
The sites that MARK never reports are now called **excluded sites**. These are known artefact positions: rCRS
310, 16182 and 16183.
- **Region annotation:** `linearized_regions.bed` and the VCF `RegionType` annotation use `Excluded_Site` (the
  homopolymer label is still `HP_Region`).
- **Sample report:** says "excluded site; not called", and its review section is headed "BELOW THE CALL THRESHOLD
  OR AT AN EXCLUDED SITE".
- **Documentation and launcher:** the README, the clean-VCF description in the launcher and the run headers use the
  same term.
- **Version:** 1.4.1 in `MARK.sh`, `MARK-I.sh`, the launcher title and `conda-recipe/meta.yaml`.

### Compatibility
A regions BED written for v1.4 or earlier, using the old label, still works. Each run's working copy of the BED is
converted to `Excluded_Site` before annotation, so those sites are still excluded and only the new label appears
in the VCFs and reports.

### Validation
- **Output comparison** against the v1.4 reference run (MARK_run_data/v1.4_20261014), on NIST C · MinION 2025 HYP
  v5.0 (3 libraries) and NIST B · Illumina 2025 (1 library): every VCF has the same records. Only the header lines
  that describe the label differ.
- **Synthetic test:** variants at rCRS 16182 and 310 are removed from the clean set with both the v1.4.1 BED and
  the old v1.4 BED.

---

<a id="v14"></a>

## MARK v1.4 release notes

v1.4 follows v1.2.10. Its one pipeline change is the variant caller: **freebayes** replaces bcftools
mpileup/call, in both MARK.sh (nanopore) and MARK-I.sh (Illumina). This was tested as v1.2.11, which was never
released. (v1.3 is a separate, parked branch with a quality-aware Illumina merge; it is not part of this release.)

### Why freebayes

bcftools judges each position on its own. freebayes first realigns the reads over the whole local stretch,
then calls.

In long C runs, Illumina read 2 puts stray T's at varying places. In the PowerSeq DNA007 sample, that left
rCRS 460T>C (in the 456–463 C run) uncalled on Illumina, although every MinION library calls it. freebayes
calls it, and nothing else changed:
- in all 63 libraries of the reference data set, the other calls are identical to v1.2.10;
- every variant position is called identically in every library of every sample.

### Changes

- **Variant calling** (both pipelines, both the review and the final BAM):

  ```
  freebayes -f ref -p 1 -E -1 --min-mapping-quality MAPQ_MIN --min-base-quality BASEQ_MIN [--limit-coverage N]
    | bcftools norm -m -any --atomize | bcftools norm -m -any
    | bcftools view -i 'GT="alt" && ALT!="*"' | bcftools sort
  ```

  - **Settings:** one copy (`-p 1`), as `--ploidy 1` before. Mapping and base quality 20, as before. No depth
    limit by default (`PILEUP_MAX_DEPTH=0`; a value above 0 becomes `--limit-coverage`).
  - **Output:** calls are split to one record per base, and only the allele the genotype calls is kept.
  - **After the raw VCF, nothing changes:** the `SAFETY_QUAL`/`MIN_DEPTH` filter, the homopolymer and
    blacklist annotation, `STRICT_QUAL` for the SNP and clean sets, and the sample report.
- **`FREEBAYES`** names the program (default `freebayes`). If it is not on `PATH`, MARK looks for a conda
  environment named `mark_fb`. If none is found, the run stops at start-up with the install command.
- **Run summary** records `VARIANT CALLER : freebayes <version> (<path>)`.
- **Launcher:** version 1.4 in the title and header. A FREEBAYES row in the parameter table, with help. The
  launcher finds freebayes by itself and logs which copy it uses. The help for `SAFETY_QUAL` and
  `STRICT_QUAL` explains the freebayes QUAL scale.

### Installing freebayes

```bash
conda create -n mark_fb -c conda-forge -c bioconda freebayes=1.3.10 -y
```

freebayes is not a dependency of the conda package. Its bioconda builds need htslib 1.21 and an older
libdeflate. On linux-64 and osx-arm64, those cannot share an environment with samtools/bcftools 1.24 and
fastp 1.3.7, the versions MARK is validated with. (They can on osx-64.) A separate environment keeps every
validated version unchanged.

### QUAL thresholds

The defaults are unchanged: `SAFETY_QUAL` 20, `STRICT_QUAL` 60. QUAL is now freebayes' score, which is not
comparable with bcftools'. It is about 20–35 per clean supporting read, so 60 is about 2–3 reads. The 10-read
`MIN_DEPTH` floor (QUAL about 200) is the stricter filter.

A calibration on 63 libraries found no confident freebayes noise call outside homopolymers:
- not at the 5% evidence floor;
- not at a 1% floor;
- not at any depth from 5x to 100x.

The confident noise that exists (16182A>C, 16183A>C, 310T>C) is already removed by the blacklist and the
homopolymer filter. Raising `STRICT_QUAL` would gain nothing and could drop real calls whose reads partly
disagree, such as DNA007 460.

### Validated with

- freebayes 1.3.10, in its own environment;
- samtools and bcftools 1.24, fastp 1.3.7, fastplong 0.7.1, cutadapt 5.2, minimap2 2.31, bwa-mem2 2.3,
  fastqc 0.12.1.

---

<a id="v1210"></a>

## MARK v1.2.10 — release notes

v1.2.10 changes one thing since **v1.2.9**: **bcftools no longer has a read-depth limit when it calls
variants.**

| | changed? |
|---|---|
| bcftools depth limit (both pipelines) | removed; see below |
| reads, trimming, merging, alignment, filters | unchanged |
| variant-calling filters | unchanged |

### 1. No depth limit in variant calling

**What was wrong.** bcftools mpileup reads at most `PILEUP_MAX_DEPTH` reads at each position. The script
default was 100,000.
- **What a pileup does at the limit:** across the whole BAM, once the limit is reached bcftools skips reads,
  and positions further along can be left with almost none.
- **The launcher made it worse:** it showed 1,000,000, but it passes a value only when its box is ticked, so
  launcher runs used 100,000 as well.

**What it did.**
- **Calls lost on a deep run:** on a re-sequencing run of DNA007 (three replicates, 3–4.7 million read pairs
  each, amplicons at 98,000–365,000 reads), two replicates lost real calls: 16129G>A, 16223C>T, 73A>G and
  263A>G. At those sites they had 1–2 reads where there were over 100,000. The nanopore data and the first
  Illumina run call all four at 99–100%.
- **Depth under-counted on an earlier run:** on the first Illumina run, the limit was reached without losing
  any call. The calls and their QUAL were the same, but the read depth behind 1–5 calls per library was
  under-counted by up to 10% (for example 2800M 152T>C: 98,802 reads, against 110,104 with no limit).

**What v1.2.10 does.**
- **No limit by default:** `PILEUP_MAX_DEPTH` now defaults to **0, meaning no limit**, in both scripts and in
  the launcher.
- **How "no limit" is passed:** the scripts give bcftools its own maximum (2,147,483,647) for both `-d` and
  `--max-idepth`. bcftools' `-d 0` also lifts the limit, but `--max-idepth 0` switches indel candidates off,
  so 0 is never passed to bcftools.
- **New log line:** every library's log gives the highest read depth in its final BAM (`[Depth] ...`).
- **New warning:** if a limit is set and any position in the final BAM reaches it, the library is listed in
  `WARNINGS.txt`, check `depth at or above PILEUP_MAX_DEPTH`, and the log says to rerun above that depth.
- **Launcher:** the `PILEUP_MAX_DEPTH` default is 0, and its help card explains the setting.

**Memory and time.** No limit costs little on these panels. The DNA007 re-sequencing run, at up to 365,000
reads per position, ran normally.

### 2. Validation

All runs used the launcher's auto cycle with every parameter at its default.

| run | final BAMs | calls |
|---|---|---|
| DNA007 re-sequencing run, 3 replicates | — | all 11 DNA007 calls in every replicate, at 99.3–100%, the same set as the first Illumina run and three nanopore replicates (16129G>A, 16223C>T, 16391G>A, 16519T>C, 73A>G, 152T>C, 199T>C, 204T>C, 207G>A, 250T>C, 263A>G) |
| first Illumina run, 4 libraries | identical to v1.2.9 | the same calls with the same QUAL; read depth behind 1–5 calls per library higher by up to 10%; NIST B identical in every respect |
| MinION 2025 SUP v5.2, 9 libraries and the unclassified bin | identical to v1.2.9, read for read | identical to v1.2.9 |
| one DNA007 replicate with `PILEUP_MAX_DEPTH=100000` | — | the new warning fires (`312067x` against `100000x`) |

### 3. Known limitations (unchanged from v1.2.9)

- **Illumina Amp10 depends on read 2 getting past the rCRS 433–438 C-tract.**
  - **DNA007:** read 2 stops at rCRS 431 in both runs (two flowcells, two library preps). Only about 0.4% of
    its Amp10 products reach the final BAM.
  - **The deeper re-sequencing run** still gives 121–188× there, enough to pass the depth check.
  - **460T>C:** on Illumina, the few reads that reach rCRS 460 are split between C and T, so it is not
    called. Nanopore calls it at 93–100%.
- **Reported calls are SNVs.** Length variants are listed in the sample report for review.

### 4. Upgrading

- **Nothing to change.** Leave `PILEUP_MAX_DEPTH` at 0.
- **If you set a limit for memory reasons,** check `WARNINGS.txt`: a library listed there may be missing
  calls.
- **Depth figures will differ:** expect higher read depth (DP, AD) than v1.2.9 at the deepest sites of deep
  libraries.

---

<a id="v129"></a>

## MARK v1.2.9 — release notes

Everything here is new since **v1.2.6**, the previous release. It covers v1.2.7, v1.2.8 and v1.2.9,
which were developed and validated together and not released separately.

**The variant-calling settings and filters are unchanged.** v1.2.9 changes three things:
- which Illumina reads are kept;
- which aligned reads reach the final BAM;
- what a run tells you when something needs attention.

The variant calls were identical to v1.2.6 in all 61 validation libraries.

---

### 1. Illumina keeps only read pairs that merge

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

### 2. Reads that align almost nothing are removed after primer clipping

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

### 3. Warnings: an amplicon below the calling depth, or losing reads after alignment

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

### 4. The per-library sample report reports the final calls

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

### 5. Length settings: the same on both platforms where possible

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

### 6. Launcher (`MARKLaunch.py`)

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

### 7. New and changed environment variables

| variable | default | pipeline | what it does |
|---|---|---|---|
| `MIN_ALIGNED_POST` | `30` | both | template bases a read must align after primer clipping (section 2) |
| `ALIGNED_WARN_PCT` | `20` | both | flag an amplicon losing this % of its fragments to that check |
| `ALIGNED_WARN_MIN` | `10` | both | ... and at least this many fragments |
| `INSERTS_BED` | `CRM_Nested_inserts.bed` | both | amplicon inserts, for the clipping step and the depth check |
| `MIN_LEN_POST` | `90` (Illumina was `30`) | both | lower end of the product-length window |
| `MAX_LEN` | `1500` | **nanopore only** (removed from Illumina) | raw read length limit |

All other variables keep their v1.2.6 names and defaults.

### 8. Validation

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

### 9. Known limitations

- **Illumina depth over Amp10 depends on read 2 surviving the rCRS 431–438 C-tract.** A library whose
  read 2 fails there loses Amp10. The depth warning names it, and the nanopore data covers the region.
- **No mate rescue on either Illumina path.** Both reads must pass and merge.
- **Reported calls are SNVs.** Indels and length variants are listed in the sample report for review,
  not called, as in earlier releases.

### 10. Upgrading

- **Calling:** nothing to change; the variant filters are the same.
- **Illumina:** expect lower Amp10 depth, and check `WARNINGS.txt` for any amplicon below the calling
  depth.
- **Scripts that read MARK outputs:** update any that read `WARNINGS.txt` (new columns),
  `final_sample_report.tsv` (new sections and depth column), the corrected VCFs (new contig name), or
  that expect intermediate BAMs or every VCF directly in `Final_Pipeline_Results/`.

---

<a id="v126"></a>

## MARK v1.2.6 — release notes

Everything here is new since **v1.2.1**, the previous release. It covers v1.2.2, v1.2.3, v1.2.4,
v1.2.5 and v1.2.6. Per-version detail is in `CHANGELOG_v1.2.3.md`, `CHANGELOG_v1.2.4.md`,
`CHANGELOG_v1.2.5.md` and `CHANGELOG_v1.2.6.md`.

**The variant-calling settings and filters are unchanged.** What changed is how reads are trimmed,
which reads are kept, and what the run reports.

---

### 1. The headline: exact amplicon trimming is now the default

`TRIM_MODE=exact` cuts every read at its own amplicon primers, and keeps a read only if it carries the
forward **and** reverse primer of the **same** amplicon with nothing but library adapter outside them.
PCR by-product — primer dimers, read-through products and amplicon-overlap products — no longer
reaches the BAM.

```
TRIM_MODE=exact      cut every read at its own amplicon primers   (default from v1.2.6)
TRIM_MODE=standard   trim with the adapter lists                  (the pre-v1.2.4 path)
```

**What to expect when you upgrade without setting anything:**

- Read counts and depth drop, because by-product is removed. The drop is larger on nanopore than on
  Illumina, and larger in libraries carrying more by-product.
- The mapped span of a full-length product equals the product size in `CRM_Nested_products.bed`.
- Read names carry the amplicon their primer pair identified (`..._Amp7`).
- On Illumina, a pair is kept only if **both** mates pass filtering (no mate rescue). A library whose
  read 1 is a dark cluster over one amplicon can lose that amplicon; check the run summary.

**To keep the old behaviour:** `TRIM_MODE=standard`, or pick it in the launcher.

Validation: across the full validation set on both platforms, reported calls were unchanged from
v1.2.1, apart from libraries that had already failed on depth. Amplicon sizes match the BED for all
10 amplicons, and nearly every read spans its amplicon exactly.

### 2. New and changed environment variables

| variable | default | pipeline | what it does |
|---|---|---|---|
| `TRIM_MODE` | `exact` | both | `exact` or `standard` (see above). The v1.2.4–v1.2.5 names `primers` and `adapters` still work and print a note. |
| `PRIMERS_BED` | `CRM_Nested_primers_empirical.bed` | both | primer coordinates; the cutadapt primer files are written from it at run time |
| `PRIMER_ERR` | `0.15` nanopore, `0.10` Illumina | both | error rate allowed when matching a primer |
| `PRIMER_MIN_OVERLAP` | `18` nanopore, `3` Illumina | both | minimum primer match length |
| `PRIMER_WARN_PCT` | `10` nanopore, `25` Illumina | both | below this share of reads carrying a primer pair, the library is flagged (the run never stops) |
| `PRIMER_DIAG_READS` | `20000` | both | how many rejected reads are classified and kept for review, per step |
| `FASTPLONG_ADAPTER_TRIM` | `0` | nanopore | `1` re-enables fastplong's own adapter auto-trimming in exact mode |
| `MINIMAP_EXTRA` | `--end-bonus 15` in exact mode | nanopore | extra minimap2 arguments |
| `BWA_EXTRA` | empty | Illumina | extra bwa-mem2 arguments |

Removed: `PRIMER_MIN_PCT` (v1.2.4 only). It aborted a run when too few reads carried a primer pair.
v1.2.5 replaced it with a warning, so a low share never stops a run.

All pre-existing variables (`MIN_LEN`, `MIN_LEN_POST`, `MAX_LEN_POST`, `QS_MIN`, `READQ`,
`UNQUAL_PCT`, `MIN_DEPTH`, `SAFETY_QUAL`, `STRICT_QUAL`, `PILEUP_MAX_DEPTH`, `EXTRA_TRIM`,
`CUTADAPT_STALL_SECS`, `MAX_LEN`, `threads`, `ref`, `regions_bed`, `RUN_NAME`, `OUTPUT_DIR`) keep
their v1.2.1 names and defaults.

### 3. Command-line changes, step by step

#### Illumina (`MARK-I.sh`)

| step | v1.2.1 | v1.2.6 |
|---|---|---|
| merge | `fastp --merge` keeping merged reads, unmerged pairs and rescued single mates | same call, but in exact mode `--unpaired1/--unpaired2` are not used: both mates must pass |
| unmerged pairs | trimmed with `cutadapt -a/-A`, `-g/-G` and aligned as pairs | v1.2.3: not aligned at all. v1.2.4+ in exact mode: kept only when read 1 and read 2 start with the two primers of one amplicon (`cutadapt --pair-adapters --action=none`) |
| second merge | – | leftover pairs are 3'-quality-trimmed and merged again: `fastp --merge --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --disable_quality_filtering --length_required 30` |
| merged reads | `cutadapt -a file:<adapters>` then `-g file:<adapters>` | exact mode: `cutadapt -g file:<anchored primers> --revcomp --action=retain --rename '{id}_{adapter_name}'`, entries `^FORWARD...REVERSE$` |
| alignment | `bwa-mem2 mem` | same; in exact mode merged reads and primer-verified pairs are aligned separately and combined with `samtools merge` |

#### Nanopore (`MARK.sh`)

| step | v1.2.1 | v1.2.6 |
|---|---|---|
| length safeguard | `fastplong --length_limit 1500` | same, plus `--disable_adapter_trimming` in exact mode (fastplong otherwise guesses an adapter and trims 10 bases past it, cutting into the primer; this also switches off its chimera splitting) |
| adapter / primer trimming | `cutadapt -b file:<adapters> --times 2` | exact mode: `cutadapt -g file:<linked primers> --revcomp --action=retain --rename '{id}_{adapter_name}' -e 0.15 -O 18`, then a per-read check that what lies outside each primer is library adapter and that no adapter sits inside the product |
| alignment | `minimap2 -ax map-ont` | same, plus `--end-bonus 15` in exact mode |

#### Both

- v1.2.2: tool steps run in a plain shell (`bash -c`) instead of a login shell (`bash -lc`). A login
  shell re-read `/etc/profile`, which on Debian and Ubuntu reset `PATH` and hid the conda environment's
  `bwa-mem2`, `samtools` and `bcftools`.

### 4. New outputs

- `WARNINGS.txt` in the run folder: one line per flagged library (its primer-pair share, the flag
  level, reads kept, and the review file). A banner also appears on screen, in the library log, and at
  the end of the run.
- `<library>_primer_rejects_review*.fastq`: up to `PRIMER_DIAG_READS` rejected reads per step, each
  labelled in its header with why it was rejected (`mark_reject=...`). Written for every library.
  Nanopore: one file. Illumina: `_merged`, `_R1` and `_R2`.
- `<library>_amplicon_trim.tsv`: reads kept per amplicon, and reads rejected per reason.
- Run summary: a `2c_PrimerPairShare` row per library, which becomes
  `2c_PrimerPairShare_LOW_see_WARNINGS` when flagged; the header records `TRIM_MODE` and whether
  fastplong's adapter auto-trimming was on.
- Read names in the BAM carry the amplicon (`..._Amp7`).

### 5. Launcher (`MARKLaunch.py`)

- `TRIM_MODE` is the first row of the parameter table, as an adapters/primers drop-down
  (`standard` / `exact`). Before this, the mode could not be chosen from the launcher.
- `FASTPLONG_ADAPTER_TRIM` is listed too (greyed out when the Illumina script is selected).
- Every parameter has a **"?"** beside its name. Clicking it explains what the setting does, which
  pipeline reads it, and its default in the selected script.

### 6. Tool versions

The conda recipe pins minimum versions to the set this release was validated with: fastp 1.3.7,
fastplong 0.7.1, cutadapt 5.2, bwa-mem2 2.3, minimap2 2.31, samtools 1.24, bcftools 1.24,
FastQC 0.12.1.

On macOS, cutadapt in multi-core mode can occasionally stop writing output part way through a large
nanopore file. The watchdog (`CUTADAPT_STALL_SECS`, default 600 s) catches this and reruns the step
single-core; the output is identical.

### 7. Known limitations

- **Depth is lower in exact mode**, by design: by-product is removed. Expect a noticeable drop in
  read count and depth, larger on nanopore than on Illumina.
- **No mate rescue in exact mode.** If read 1 is a dark cluster over one amplicon, that amplicon can
  lose almost all of its depth even though read 2 is fine. Check the per-amplicon counts in the run
  summary, or use `TRIM_MODE=standard` for such a library.
- **fastplong's own quality filter** (its Q15 default, which MARK does not set) removes a small share
  of whole products on good nanopore libraries, and a much larger share on low-quality ones.
- Reported calls are SNVs; indels are not reported in the main VCFs, unchanged from earlier releases.

### 8. Upgrading

- Nothing to change for calling: the variant filters are the same.
- Expect fewer reads and lower depth than v1.2.1 on the same data, because by-product is removed,
  alongside cleaner amplicon boundaries and a higher share of reads that are a whole product.
- To reproduce v1.2.1-style read retention, set `TRIM_MODE=standard`.

## MARK v1.2.6 — clearer trimming-mode names; fastplong adapter auto-trimming off in exact mode

This release includes everything in v1.2.5 (see `CHANGELOG_v1.2.5.md`). Two things change what a run
does by default: exact amplicon trimming is now the default mode, and fastplong's own adapter
auto-trimming is switched off in that mode. Variant-calling settings and filters are unchanged, and
`TRIM_MODE=standard` restores the pre-v1.2.4 behaviour.

### Exact amplicon trimming is now the default

`TRIM_MODE` defaults to `exact` from this release (it was `adapters`/`standard` in v1.2.4 and v1.2.5).
Exact trimming was validated against the v1.2.1 release across the validation set with identical calls, and it
is what removes PCR by-product from the BAM. To keep the pre-v1.2.4 behaviour, set
`TRIM_MODE=standard`, or choose it in the launcher.

What changes for a run that does not set `TRIM_MODE`:

- Reads are cut exactly at their own amplicon primers, and a read is kept only if it carries both
  primers of one amplicon. PCR by-product (primer dimers, read-through and amplicon-overlap products)
  no longer reaches the BAM, so **read counts and depth drop**, more so on nanopore than on Illumina.
- Illumina pairs are kept only when both mates pass filtering (no mate rescue). A library whose read 1
  is a dark cluster over one amplicon can therefore lose that amplicon; check the run summary.
- Read names carry the amplicon their primer pair identified (`..._Amp7`).

### Trimming modes renamed

Both modes remove the library adapters, and in both the primer footprints are soft-clipped before
calling, so "adapters" versus "primers" read like a choice between two things every run needs. The
modes differ in how a read is cut and which reads are kept, and are now named for that:

| v1.2.4–v1.2.5 | v1.2.6 | what it does |
|---|---|---|
| `TRIM_MODE=adapters` | `TRIM_MODE=standard` | adapter sequences from the adapter list are cut wherever found; any read long enough afterwards is kept |
| `TRIM_MODE=primers` | `TRIM_MODE=exact` (**now the default**) | exact amplicon trimming: a read is kept only with both primers of one amplicon, and cut exactly at them |

The old names still work and print a one-line note. Results are unchanged. The run summary header
now records the mode with a short description, for Illumina as well as nanopore.

### fastplong adapter auto-trimming off in exact mode

By default fastplong guesses a library adapter for each file and trims it, plus 10 more bases. In
exact mode those 10 bases are the start of the primer, so the read then fails the primer check or
is kept clipped short. The primer step already removes everything outside the primers, so in exact
mode the auto-trimming adds nothing and costs good reads.

| | v1.2.5 | v1.2.6 |
|---|---|---|
| exact mode | fastplong adapter auto-trimming on | **off** (`--disable_adapter_trimming`) |
| standard mode | on | on (unchanged) |
| new setting | – | `FASTPLONG_ADAPTER_TRIM` (default `0`); `1` turns the auto-trimming back on in exact mode |

- The setting is also in the launcher's parameter table, as `FASTPLONG_ADAPTER_TRIM`.
- The launcher can now choose the trimming mode. `TRIM_MODE` is the first row of its parameter table,
  as a drop-down (standard / exact) that follows the script default until overridden. Before this,
  runs started from the launcher could not select the mode at all.
- The run summary header records `TRIM_MODE` and whether fastplong auto-trimming was on.
- fastplong's chimera splitting is not set by MARK. Note that `--disable_adapter_trimming` turns it
  off as well, because splitting looks for adapter inside a read: in exact mode no read is split.
  This is harmless here - a chimera carries adapter inside the product, which the adapter-junction
  check rejects - and switching splitting off on its own made no difference to real libraries.
- `--trimming_extension 0` was tested as a way to keep the auto-trimming without the 10-base
  overshoot. It does not help: the reads kept either way end up spanning their amplicon less exactly,
  and the reads it adds are markedly less exact than the ones it drops. fastplong's adapter detection
  is not precise enough on noisy reads to stop cleanly at the adapter.

### Launcher: help for every parameter

Each parameter in the launcher's table now has a "?" beside its name. Clicking it opens a short card
that says what the setting does, which pipeline uses it (nanopore, Illumina or both), and its default
in the selected script. The card closes with Close, Escape, or a click elsewhere.
`DISCARD_WARN_PCT` is defined in both scripts but read by no step; its card says so.

### Evidence

The auto-trimming only acts where fastplong's guess finds an adapter: in the validation data, the 9
libraries basecalled with the fast model, and in unclassified barcode bins; it touched no reads in
libraries basecalled with the high-accuracy model. With it off in exact mode, on the affected
libraries:

- reads kept roughly doubled, and the lowest-covered amplicon gained the most
- reads that did not span one whole product fell to a handful per library
- calls that had been lost against v1.2.1 in a low-yield library came back
- the error rate rose only marginally, and no other calls changed

Details: `development/test_builds/MARK_FastplongDefaults_Test_20260922`.

### Verified

- **Names:** `TRIM_MODE=exact` and the old `TRIM_MODE=primers` give identical results on the same
  library (8,327 reads); `standard` on that library keeps 13,042. `TRIM_MODE=foo` stops with
  "TRIM_MODE must be standard or exact".
- **fastplong auto-trimming:**

On one low-quality library and one good nanopore library:

| case | low-quality library, reads | good library, reads |
|---|---|---|
| exact mode, default | 4,107 | 8,327 (unchanged) |
| exact mode, `FASTPLONG_ADAPTER_TRIM=1` | 2,256 (identical to v1.2.5) | 8,327 |
| standard mode | 7,812 (identical to v1.2.5) | 13,042 |

### Also

- Version 1.2.6 in `MARK.sh`, `MARK-I.sh`, `MARKLaunch.py` and `conda-recipe/meta.yaml`.

---

<a id="v125"></a>

## MARK v1.2.5 — primer-pair check warns and never stops the run

This release includes everything in v1.2.4 (see `CHANGELOG_v1.2.4.md`). Nothing changes in trimming,
alignment, calling or filtering: the same reads reach the BAM as in v1.2.4. v1.2.4 remains available
as a fallback.

### What changed (TRIM_MODE=exact only)

In exact mode, reads without a primer pair of one amplicon are removed, exactly as before. What
changes is what happens when a library has only a small share of such reads.

| | v1.2.4 | v1.2.5 |
|---|---|---|
| reads without a primer pair | removed | removed (unchanged) |
| share below `PRIMER_MIN_PCT` (20%) | the whole run stops | setting removed: the run never stops |
| share below `PRIMER_WARN_PCT` | one line in the log | flagged prominently (below) |
| `PRIMER_WARN_PCT` default | 25% nanopore, 40% Illumina | **10% nanopore, 25% Illumina** |
| rejected reads | classified on a sample, then deleted | the same sample is kept, labelled, in a review FASTQ for **every** library |

Why: the share of reads with a primer pair varies widely between good libraries, because most of
what is sequenced is PCR by-product. Across the validation set that share ranged from a few per cent
to more than half, so the old 20% stop would have halted most nanopore runs. A stop also ended the
whole run, not just the library. The kept reads of a low-share library are clean, primer-verified
products, so there is no reason to discard them.

#### How a flagged library is shown

- A banner in the terminal and in the library's log.
- `WARNINGS.txt` at the top of the run folder, one line per flagged library: its share, the flag
  level, the reads kept and the review file.
- The run summary row `2c_PrimerPairShare` becomes `2c_PrimerPairShare_LOW_see_WARNINGS`.
- A second banner at the end of the run lists every flagged library.

#### Review files (every library)

Up to `PRIMER_DIAG_READS` (default 20,000) rejected reads per rejecting step, as they entered the
primer step, with the reason in the read header (`mark_reject=...`):

- Nanopore: `<library>_primer_rejects_review.fastq`, holding reads rejected by
  - the primer step: no primer found, only one primer of one amplicon, or primers from two amplicons
  - the adapter check: non-adapter sequence outside the primers, or adapter inside the product
- Illumina:
  - `<library>_primer_rejects_review_merged.fastq`: merged reads that are not one whole product
  - `<library>_primer_rejects_review_R1.fastq` and `_R2.fastq`: unmerged pairs whose mates do not
    start with the two primers of one amplicon

At the defaults each file is a few megabytes. Illumina's pair check now writes rejected pairs to
a file instead of discarding them (`--untrimmed-output` in place of `--discard-untrimmed`). The
pairs that are kept are the same.

### Fixed

- **Reject reasons in the run summary.** In v1.2.4, rejected reads with no primer match at all were
  never counted, so "no primer found" was missing from `<library>_amplicon_trim.tsv`. Those reads
  were spread across "only one primer" and "primers from two amplicons". They are now counted
  (a few per cent of the rejected reads in a low-yield library). Which reads are kept is unaffected.

### Verified

On nanopore and Illumina libraries, including a low-yield one, v1.2.5 keeps exactly the same reads as
v1.2.4 (identical read names in the final BAM), and the snps and clean VCFs are identical.

- A library whose primer-pair share fell below the flag level was flagged and still processed.
- A library above the level was flagged only when the level was raised for the test.

### Also

- Version 1.2.5 in `MARK.sh`, `MARK-I.sh`, `MARKLaunch.py` and `conda-recipe/meta.yaml`.

---

<a id="v124"></a>

## MARK v1.2.4 — exact amplicon trimming

This release also includes the changes made in v1.2.2 and v1.2.3 (see `CHANGELOG_v1.2.3.md`). The
variant-calling settings and filters are unchanged.

### A second trimming mode: `TRIM_MODE`

```
TRIM_MODE=standard   trim with the adapter lists (the pre-v1.2.4 path)
TRIM_MODE=exact      cut every read at its own amplicon primers
```

In v1.2.4 and v1.2.5 these values were named `adapters` and `primers`, and `adapters` was the default.
v1.2.6 renamed them (the old names still work) and made `exact` the default.

Trimming with adapter lists can never be exact by construction. An adapter entry that stops one base
short of the junction leaves a base on the read; an entry that contains primer bases cuts into the
amplicon. In exact mode a read is kept only if it carries the forward **and** reverse primer of the
**same** amplicon, with nothing but library adapter outside them, and is cut exactly at the outer
primer ends. The mapped span of a full-length product is then the product size in
`CRM_Nested_products.bed` by construction.

- The cutadapt primer files are written at run time from `CRM_Nested_primers_empirical.bed` and the
  reference, so they cannot drift from the coordinates the rest of the pipeline uses.
- The primers stay on the read (`--action=retain`) on purpose: the product coordinates include them,
  and the existing primer-footprint step soft-clips them in the BAM, so no primer base reaches the
  caller.
- Each read name gains the amplicon its primer pair identified (`..._Amp7`), so amplicon assignment no
  longer depends on alignment coordinates alone.
- The run summary reports, per amplicon, reads kept and reads rejected with the reason (no primer pair,
  only one primer, primers from two amplicons, non-adapter sequence outside the primers, adapter inside
  the product).

#### Illumina (`MARK-I.sh`), exact mode

1. fastp merge as before, but **no mate rescue**: a fragment is kept only if both mates survive.
2. Pairs that did not merge are quality-trimmed at the 3' end (`--cut_right` Q20) and merged again.
   (`--cut_right` in the first pass is worse on every amplicon, so it is used only here.)
3. Merged reads must be one whole product: `cutadapt -g file:<anchored> --revcomp --action=retain`,
   where each entry is `^FORWARD...REVERSE$`. Anchoring is what rejects read-through products; without
   it a product spanning two amplicons would be trimmed down into a perfect-looking shorter one.
4. Pairs that still did not merge are kept only if read 1 and read 2 start with the two primers of one
   amplicon (`--pair-adapters`), and are aligned as a pair.
5. The `-a`/`-g` adapter-list passes are not used in this mode.

#### Nanopore (`MARK.sh`), exact mode

The `-b` adapter-list pass is replaced by the linked-primer pass (`-e 0.15 -O 18`), followed by a filter
that requires the sequence outside each primer to be the library adapter junction, and rejects a product
that contains adapter **inside** it (a chimera whose junction falls between two primer matches).
`minimap2` gains `--end-bonus 15`: with exact trimming the read ends are real template, and without it
the aligner soft-clips error-containing primer ends. The flag markedly raises the share of reads that
span their amplicon exactly, and leaves calls unchanged.

#### Guard

If too few reads carry a matching primer pair, the run flags the library, so a kit or primer change is
visible instead of silently discarding everything. Most of what is sequenced is by-product, so the
share is a minority of reads by nature, and varies widely between libraries. In v1.2.4 a low share
aborted the run; v1.2.5 replaced that with a warning (see `CHANGELOG_v1.2.5.md`).

### Validation

The full validation set, both platforms, against the validated v1.2.1 runs:

- Illumina: reported calls identical, with allele fractions essentially unchanged.
- Nanopore: calls identical except in libraries that had already failed on depth.
- Amplicon sizes match the BED for all 10 amplicons on both platforms, and nearly every read spans its
  amplicon exactly.
- No library adapter survives trimming, and no read-through product reaches the BAM.
- Cost: read count and mean depth fall, more so on nanopore, and what is dropped is by-product. Where
  a minimum depth falls, most of the lost depth was read-through product rather than genuine template.

---

<a id="v123"></a>

## MARK v1.2.3 — changes since v1.2.1

This release also includes the changes made in v1.2.2. The variant-calling settings and filters are unchanged.

### Illumina: merged reads only (`MARK-I.sh`)

```
v1.2.1   merged reads (single-end)  +  unmerged pairs trimmed with cutadapt -a/-A, -g/-G and aligned as pairs
v1.2.3   merged reads (single-end) only
```

v1.2.0 began aligning the pairs fastp cannot merge. Part of that set is longer PCR product that runs from the forward primer of one amplicon to the reverse primer of a later one. At 300–550 bp it is too long to merge, and it appears in IGV as large-insert (red) pairs. In v1.2.3 only merged reads are aligned, so every aligned read is one collapsed fragment.

- fastp still writes unmerged pairs to `<sample>_unmerged_R1.fastq` / `_R2.fastq`, so they can be counted and inspected. The log reports how many there were.
- Reads whose mate failed quality filtering (`--unpaired1/--unpaired2`) are still kept as single-end reads, as in v1.2.1.
- The run summary no longer has a `2b_TrimPairs` row. `1_MergeFilter(Reads)` now counts merged reads only.
- Depth is lower where most fragments do not merge, typically across length-variable homopolymers. Check per-position depth there before relying on a call.

### Tool versions

The pipeline has been validated with:

| tool | version |
|---|---|
| fastp | 1.3.7 |
| fastplong | 0.7.1 |
| cutadapt | 5.2 |
| bwa-mem2 | 2.3 |
| minimap2 | 2.31 |
| samtools | 1.24 |
| bcftools | 1.24 |
| FastQC | 0.12.1 |

The conda recipe now sets these as minimum versions.

On macOS, cutadapt in multi-core mode can occasionally stop writing output part way through a large nanopore file, whatever the cutadapt version. The existing watchdog (`CUTADAPT_STALL_SECS`, default 600 s) catches this and reruns the step single-core. The output is identical; the step just takes longer.

### Linux: tool steps run in a plain shell (v1.2.2, `MARK.sh`, `MARK-I.sh`)

```
v1.2.1   run_log bash -lc "{ ... }"
v1.2.2   run_log bash -c  "{ ... }"
```

A login shell re-reads `/etc/profile`. On Debian and Ubuntu that file sets `PATH` outright, which drops the activated conda environment, so the alignment and calling steps could not find `bwa-mem2`, `samtools` or `bcftools`. On macOS it moved `/usr/local/bin` ahead of the environment, so a system copy of a tool could be used instead of the environment's. The pipelines no longer start login shells.

---

<a id="v121"></a>

## MARK v1.2.1 — changes since v1.1.8

This release also includes the changes made in v1.1.9. The variant-calling settings and filters are unchanged. What changed is how reads are filtered, trimmed and aligned before calling: more real template reaches the caller, and adapter and primer sequence is no longer counted as sample.

### Adapter lists (both pipelines)

| | v1.1.8 | v1.2.1 |
|---|---|---|
| Illumina (`MARK-I.sh`) | `Updated_Adapter_Primer_List_Cutadapt_cleaned.txt`, 744 entries, shared | `MARK_Adapter_List_Illumina.txt`, 200 entries |
| Nanopore (`MARK.sh`) | the same shared list | `MARK_Adapter_List_ONT.txt`, 207 entries (Illumina + Nanopore adapters) |

- **Mitochondrial sequence removed.** The old list held 320 entries of mtDNA control-region sequence, plus several adapter–mtDNA chimeras. cutadapt `-a` removes a match *and everything after it*, and `-g` removes a match *and everything before it*. So each mtDNA entry cut every read containing it at a fixed position in the genome, and the remainder often failed the length filter. No entry now has a ≥15 bp exact match to the reference. `MARK_Adapter_List_v2.rebuild_log.txt` records every edit.
- **Index barcodes removed.** These 8 bp sample indices are read in a separate index cycle and never occur in the insert. With `--overlap 5`, a chance 5 bp match could trim a real read end.
- **Nanopore chemistry removed from the Illumina list.** Illumina libraries cannot contain it.
- The dashboard now switches to the matching list when you change pipelines, unless you have chosen your own.

### Nanopore adapter trimming (`MARK.sh`)

```
v1.1.8   cutadapt -a file:LIST ...   then   cutadapt -g file:LIST ...
v1.2.1   cutadapt -b file:LIST --times 2 ...
```

These amplicons carry adapter at both ends of the read. `-a` treats every match as a 3′ adapter, so an adapter matched at the start of a read removed the entire read, which `--minimum-length` then discarded. `-b` treats a match that includes the read's first base as a 5′ adapter, keeping the sequence after it, and any other match as a 3′ adapter. `--times 2` lets one read lose an adapter from each end. Length floors are unchanged (`MIN_LEN=90`, `MIN_LEN_POST=90`).

### Illumina read handling (`MARK-I.sh`)

```
v1.1.8   fastp --merge --include_unmerged --merged_out merged.fq --length_required 90 ...
v1.2.1   fastp --merge --trim_poly_g --merged_out merged.fq \
               --out1 unmerged_R1.fq --out2 unmerged_R2.fq \
               --unpaired1 rescued_R1.fq --unpaired2 rescued_R2.fq --length_required 90 ...
```

- **`--trim_poly_g`** removes poly-G tails, which two-colour Illumina chemistry produces when a cluster gives no signal. Left in, they align with most of the read soft-clipped and add depth that is not sample.
- **Unmerged pairs stay paired.** `--include_unmerged` placed them in the merged file as two unrelated single-end reads, so where mates overlapped, one fragment was counted twice. Now each unmerged pair goes through four steps:
  1. it is written to `--out1`/`--out2`;
  2. it is trimmed in paired mode (`cutadapt -a/-A`, then `-g/-G`);
  3. it is aligned as a pair (`bwa-mem2 mem ref R1 R2`);
  4. it is combined with the merged reads (`samtools merge`).

  After primer trimming, the BAM passes through `samtools sort -n | samtools fixmate | samtools sort`, so mate positions stay correct and the pileup counts each fragment once.
- **A good mate is no longer discarded with a bad one.** By default fastp drops the whole pair when either read fails filtering. `--unpaired1`/`--unpaired2` keep the read that passed, and it continues as single-end. This matters when a fault affects only one read of the pair, such as poly-G in R1.
- **The 90 bp floor applies once, to the raw fragment.** v1.1.8 reapplied it at every cutadapt step and at the post-trim filter, which discarded valid fragments once their adapter and primer were removed. After trimming, the floor is now `MIN_LEN_POST=30`, a minimum for reliable mapping (default was 90).

### Primer trimming (both pipelines)

**v1.1.8:** each read was assigned to the amplicon it overlapped most, then cut to fixed, non-overlapping tiles.

**v1.2.1:** each read end is compared with the known PCR product ends, within `PRIMER_TOL = 10` bp.
- An end at a product start has that product's forward primer removed.
- An end at a product end has its reverse primer removed.
- Everything between the primers is kept, whichever amplicons it spans.
- A read with neither end at a known product end falls back to the old rule and is clipped to a single insert.

This gives four improvements:

- **Fragments spanning two amplicons keep all their template.** Examples are a forward primer paired with the neighbouring amplicon's reverse primer. Previously, such a fragment was clipped to one amplicon.
- **Primer is removed by position, not by amplicon assignment.** Primer sequence always matches the reference, so counting it as coverage can hide a real variant under its footprint.
- **The insert table matches the ten inserts published for the PowerSeq CRM Nested System.** Adjacent inserts therefore overlap as they do in the kit, and a base covered by two amplicons counts reads from both.
- **Reads outside every amplicon are dropped as off-target.** Previously they crashed the trimming step.

The tolerance matters because the number of primer bases left at the start of a read after adapter trimming varies, and differs between Illumina and nanopore. The closest termini of two different amplicons are 13 bp apart, so a tolerance of 10 cannot match the wrong amplicon.

#### Amplicon coordinates

The inserts are those published by Vinueza-Espinosa et al. 2023 (*Electrophoresis* 44:1423-1434). Promega does not publish the primer sequences, so the product spans are derived from the fixed read termini seen when sequencing this kit. The pipelines run against a linearized reference, `NC_012920.1_linearized`, where `rCRS = (lin + 8284) mod 16569`. Coordinates below are 1-based and inclusive.

| Amplicon | Insert (rCRS) | Insert (linearized) | PCR product (linearized) |
|---|---|---|---|
| Amp1 | 16013–16126 | 7729–7842 | 7702–7868 |
| Amp2 | 16116–16225 | 7832–7941 | 7810–7964 |
| Amp3 | 16223–16408 | 7939–8124 | 7913–8149 |
| Amp4 | 16387–16486 | 8103–8202 | 8079–8225 |
| Amp5 | 16474–30 | 8190–8315 | 8166–8337 |
| Amp6 | 16555–152 | 8271–8437 | 8249–8465 |
| Amp7 | 136–257 | 8421–8542 | 8394–8578 |
| Amp8 | 246–364 | 8531–8649 | 8503–8674 |
| Amp9 | 342–436 | 8627–8721 | 8602–8745 |
| Amp10 | 429–592 | 8714–8877 | 8687–8904 |

The same coordinates ship as BED files (0-based, half-open): `CRM_Nested_inserts.bed`, `CRM_Nested_products.bed` and `CRM_Nested_primers_empirical.bed`. The pipelines carry their own copy internally and do not read these files.

For Amp6 and Amp9, adapter trimming usually removes the first few primer bases, so reads start at 8258 and 8608 rather than at the product starts above. The trimmer's internal table records those positions, and `PRIMER_TOL` covers both cases.

#### How the trimmer handles each kind of read

Linearized coordinates, 1-based.

| Read | What it is | Kept |
|---|---|---|
| 8258–8465 | Amp6 product, primer partly removed by adapter trimming | 8271–8437 |
| 8249–8465 | Amp6 product, full primer still attached | 8271–8437 |
| 8079–8337 | Amp4 forward primer to Amp5 reverse primer (hybrid) | 8103–8315 |
| 8249–8337 | Amp6 forward primer to Amp5 reverse primer (hybrid) | 8271–8315 |
| 8687–8723 | Amp10 read that stops early | 8714–8723 |
| 8727–8745 | read lying entirely inside Amp9's reverse primer | dropped |

### Variant calling (both pipelines)

`bcftools mpileup` now also receives `--max-idepth $PILEUP_MAX_IDEPTH`, which defaults to `PILEUP_MAX_DEPTH`. bcftools applies a separate depth cap of 250 to indel candidates. Raising `-d` alone left that cap in place, so no indels were proposed above 250×. This affects only the `annotated_all` and `qual_filtered` review VCFs; the `snps` and `clean` outputs remain SNP-only.

### Running the pipeline

- **cutadapt can no longer hang a run.** Its multi-core mode can deadlock on macOS. Every cutadapt call now watches its output file. If the output stops growing for `CUTADAPT_STALL_SECS` seconds, the call is rerun with `--cores 1`, and a second stall stops the pipeline with an error. Output is identical either way, because cutadapt keeps reads in input order.
- **The scripts use their own data files.** The reference, regions BED and adapter list are found beside the script, following symlinks, so command-line runs work from any directory. Explicit `ref=`, `regions_bed=` and `ADAPTER_FILE=` values still take priority. The resolved paths are recorded in `run_summary.txt`.
- **The run folder is created beside the input.** Previously a command-line run put it in the current directory, and a folder browsed in the dashboard could place it inside the input. `OUTPUT_DIR` sets another location, which is created if missing. `RUN_NAME` names the folder, and the dashboard's Custom Run Name field now works. A run will not write into an existing folder that is not empty.

### New and changed settings

| Setting | Default | Pipeline |
|---|---|---|
| `CUTADAPT_STALL_SECS` | 600 | both — new |
| `PILEUP_MAX_IDEPTH` | same as `PILEUP_MAX_DEPTH` | both — new |
| `OUTPUT_DIR` | the folder holding the input | both — new |
| `RUN_NAME` | auto-generated | both — new |
| `MIN_LEN_POST` | 30 (was 90) | `MARK-I.sh` |

### Upgrading from v1.1.8

- **Custom adapter lists:** pass them with `ADAPTER_FILE=` or the dashboard's adapter field. A list placed next to the input is no longer picked up, and the old file name is no longer searched for.
- **Output location:** run folders now appear beside the input, not in the directory the script was started from.
- **`run_summary.txt`:**
  - Illumina reports `2_Trim3(merged)` and `2b_TrimPairs` separately, plus a count of rescued mates.
  - Nanopore reports a single `2_AdapterTrim` step.
  - Both record `PILEUP_IDEPTH` and the resolved reference, regions BED and adapter list paths.

Version 1.2.1 across the dashboard, both pipelines and the conda recipe.
