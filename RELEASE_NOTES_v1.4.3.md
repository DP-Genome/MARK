# MARK v1.4.3

**Primer bases are hard-clipped after alignment by default. The variant calls are identical to v1.4.2.**

## What changed
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

## Why hard clipping
- The final (`trimmed`/`cleaned`) BAMs are about 20% smaller.
- They hold only the bases that variant calling uses.
- The untrimmed baseline BAM still shows the primers for review.

## Validation
- **Hard vs soft:** on v1.4.1 with `CLIP_MODE=hard`, compared with soft clipping across 10 libraries (DNA 007 ·
  Illumina 2026 ×3, DNA 007 · MinION 2025 SUP v5.2 ×3, NIST C · MinION 2023 SUP v5.2 ×3, NIST B · Illumina 2025):
  - all 240 VCFs had identical records;
  - the trimmed BAMs had the same reads and the same aligned bases;
  - every soft clip became a hard clip.
- **v1.4.3 with default settings:** DNA 007 · MinION 2025 SUP v5.2 (3 libraries) and 2800M · Illumina 2025. All
  192 VCFs are identical to v1.4.2, the final BAMs carry hard clips only, and the run summary records CLIP_MODE=hard.
- **`CLIP_MODE=soft`:** reproduces the v1.4.2 CIGARs (e.g. 22S110M23S for a whole Amp2 product; hard: 22H110M23H).
