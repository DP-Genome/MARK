# MARK v1.6.2

**Primer bases are always hard-clipped.** The `CLIP_MODE` setting is removed: after alignment, the primer bases at
each product end are removed from the read (CIGAR `H`), which has been the default since v1.4.3. The soft-clip
option, its launcher drop-down and its help text are gone, and the run summary states the clipping as
"PRIMER CLIPPING: primer bases hard-clipped". A `CLIP_MODE` value left in the environment is ignored.

The pipelines, the calls and the heteroplasmy report are otherwise identical to v1.6.1: re-running all 25 reference
runs (63 libraries) gave the same alignments, VCFs, read counts and trimmed reads, and sample reports that differ only
in their title line ("MARK v1.6.2 sample report").
