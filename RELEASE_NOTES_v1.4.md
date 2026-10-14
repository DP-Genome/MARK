# MARK v1.4 release notes

v1.4 follows v1.2.10. Its one pipeline change is the variant caller: **freebayes** replaces bcftools
mpileup/call, in both MARK.sh (nanopore) and MARK-I.sh (Illumina). This was tested as v1.2.11, which was never
released. (v1.3 is a separate, parked branch with a quality-aware Illumina merge; it is not part of this release.)

## Why freebayes

bcftools judges each position on its own. freebayes first realigns the reads over the whole local stretch,
then calls.

In long C runs, Illumina read 2 puts stray T's at varying places. In the PowerSeq DNA007 sample, that left
rCRS 460T>C (in the 456–463 C run) uncalled on Illumina, although every MinION library calls it. freebayes
calls it, and nothing else changed:
- in all 63 libraries of the reference data set, the other calls are identical to v1.2.10;
- every variant position is called identically in every library of every sample.

## Changes

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

## Installing freebayes

```bash
conda create -n mark_fb -c conda-forge -c bioconda freebayes=1.3.10 -y
```

freebayes is not a dependency of the conda package. Its bioconda builds need htslib 1.21 and an older
libdeflate. On linux-64 and osx-arm64, those cannot share an environment with samtools/bcftools 1.24 and
fastp 1.3.7, the versions MARK is validated with. (They can on osx-64.) A separate environment keeps every
validated version unchanged.

## QUAL thresholds

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

## Validated with

- freebayes 1.3.10, in its own environment;
- samtools and bcftools 1.24, fastp 1.3.7, fastplong 0.7.1, cutadapt 5.2, minimap2 2.31, bwa-mem2 2.3,
  fastqc 0.12.1.
