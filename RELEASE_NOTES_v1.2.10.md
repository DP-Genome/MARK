# MARK v1.2.10 — release notes

v1.2.10 changes one thing since **v1.2.9**: **bcftools no longer has a read-depth limit when it calls
variants.**

| | changed? |
|---|---|
| bcftools depth limit (both pipelines) | removed; see below |
| reads, trimming, merging, alignment, filters | unchanged |
| variant-calling filters | unchanged |

## 1. No depth limit in variant calling

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

## 2. Validation

All runs used the launcher's auto cycle with every parameter at its default.

| run | final BAMs | calls |
|---|---|---|
| DNA007 re-sequencing run, 3 replicates | — | all 11 DNA007 calls in every replicate, at 99.3–100%, the same set as the first Illumina run and three nanopore replicates (16129G>A, 16223C>T, 16391G>A, 16519T>C, 73A>G, 152T>C, 199T>C, 204T>C, 207G>A, 250T>C, 263A>G) |
| first Illumina run, 4 libraries | identical to v1.2.9 | the same calls with the same QUAL; read depth behind 1–5 calls per library higher by up to 10%; NIST B identical in every respect |
| MinION 2025 SUP v5.2, 9 libraries and the unclassified bin | identical to v1.2.9, read for read | identical to v1.2.9 |
| one DNA007 replicate with `PILEUP_MAX_DEPTH=100000` | — | the new warning fires (`312067x` against `100000x`) |

## 3. Known limitations (unchanged from v1.2.9)

- **Illumina Amp10 depends on read 2 getting past the rCRS 433–438 C-tract.**
  - **DNA007:** read 2 stops at rCRS 431 in both runs (two flowcells, two library preps). Only about 0.4% of
    its Amp10 products reach the final BAM.
  - **The deeper re-sequencing run** still gives 121–188× there, enough to pass the depth check.
  - **460T>C:** on Illumina, the few reads that reach rCRS 460 are split between C and T, so it is not
    called. Nanopore calls it at 93–100%.
- **Reported calls are SNVs.** Length variants are listed in the sample report for review.

## 4. Upgrading

- **Nothing to change.** Leave `PILEUP_MAX_DEPTH` at 0.
- **If you set a limit for memory reasons,** check `WARNINGS.txt`: a library listed there may be missing
  calls.
- **Depth figures will differ:** expect higher read depth (DP, AD) than v1.2.9 at the deepest sites of deep
  libraries.
