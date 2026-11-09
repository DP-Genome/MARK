# MARK v1.7: mixture check, heteroplasmy report fix

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

## How it works
A second person's bases at the positions where the two haplotypes differ sit **on the same reads**. Sequencing errors and DNA damage also leave minor bases, but each sits on its own read. The new script `MARK_mixture_check.py` runs on the final, primer-clipped BAM:

1. A candidate is the second most common base at a position, outside homopolymer regions and excluded sites, carried by at least 1% of reads and at least 5 reads.
2. Two candidates are **linked** when they meet all of these:
   - at least 50 reads cover both positions;
   - at least 5 of those reads carry both minor bases;
   - that count is at least 70% of the smaller minor-base count;
   - the share of covering reads that carry both is at least `MIX_MIN_AF`.
3. One linked pair is enough to flag the library. The estimated share of the second haplotype is the median over the linked pairs.

## Outputs
- **WARNINGS.txt:** "evidence of a mixture (second haplotype)", with the amplicons, the estimated share and the linked bases.
- **The sample report:** a new section "MIXTURE CHECK: A SECOND HAPLOTYPE?". It is written for every library, flagged or not.
- **`<lib>_mixture_check.tsv`:** the linked pairs and their read counts, copied to `sample_reports/`.
- **Run summary and launcher:** a run-summary line, and the console banner during the run.

## Settings
- `MIX_MIN_AF` (default 0.02): the smallest share of the second haplotype that is reported. `0` turns the check off.
- `HET_MIN_AF` (unchanged, default 0.05) can still be raised or lowered.

## Validation
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
