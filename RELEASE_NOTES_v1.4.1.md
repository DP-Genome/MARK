# MARK v1.4.1

**Terminology change only. The variant calls are identical to v1.4.**

## What changed
The sites that MARK never reports are now called **excluded sites**. These are known artefact positions: rCRS
310, 16182 and 16183.
- **Region annotation:** `linearized_regions.bed` and the VCF `RegionType` annotation use `Excluded_Site` (the
  homopolymer label is still `HP_Region`).
- **Sample report:** says "excluded site; not called", and its review section is headed "BELOW THE CALL THRESHOLD
  OR AT AN EXCLUDED SITE".
- **Documentation and launcher:** the README, the clean-VCF description in the launcher and the run headers use the
  same term.
- **Version:** 1.4.1 in `MARK.sh`, `MARK-I.sh`, the launcher title and `conda-recipe/meta.yaml`.

## Compatibility
A regions BED written for v1.4 or earlier, using the old label, still works. Each run's working copy of the BED is
converted to `Excluded_Site` before annotation, so those sites are still excluded and only the new label appears
in the VCFs and reports.

## Validation
- **Output comparison** against the v1.4 reference run (MARK_run_data/v1.4_20261014), on NIST C · MinION 2025 HYP
  v5.0 (3 libraries) and NIST B · Illumina 2025 (1 library): every VCF has the same records. Only the header lines
  that describe the label differ.
- **Synthetic test:** variants at rCRS 16182 and 310 are removed from the clean set with both the v1.4.1 BED and
  the old v1.4 BED.
