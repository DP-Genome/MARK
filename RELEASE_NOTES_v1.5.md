# MARK v1.5

**Choose the variant caller: freebayes (the default), bcftools, or both, with a per-library comparison of the two.
With the default, the calls are identical to v1.4.4.**

## What changed
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

## Thresholds
Both callers go through the same filters (depth ≥ `MIN_DEPTH`, QUAL > `SAFETY_QUAL` for review, SNVs with QUAL >
`STRICT_QUAL` for calls). These were calibrated for freebayes; bcftools' QUAL is on a different scale (it tops out
around 228), so a bcftools-only or bcftools-missed variant near the threshold should be read in the alignment.

## Validation
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
