# MARK v1.2.6 — release notes

Everything here is new since **v1.2.1**, the previous release. It covers v1.2.2, v1.2.3, v1.2.4,
v1.2.5 and v1.2.6. Per-version detail is in `CHANGELOG_v1.2.3.md`, `CHANGELOG_v1.2.4.md`,
`CHANGELOG_v1.2.5.md` and `CHANGELOG_v1.2.6.md`.

**The variant-calling settings and filters are unchanged.** What changed is how reads are trimmed,
which reads are kept, and what the run reports.

---

## 1. The headline: exact amplicon trimming is now the default

`TRIM_MODE=exact` cuts every read at its own amplicon primers, and keeps a read only if it carries the
forward **and** reverse primer of the **same** amplicon with nothing but library adapter outside them.
PCR by-product — primer dimers, read-through products and amplicon-overlap products — no longer
reaches the BAM.

```
TRIM_MODE=exact      cut every read at its own amplicon primers   (default from v1.2.6)
TRIM_MODE=standard   trim with the adapter lists                  (the pre-v1.2.4 path)
```

**What to expect when you upgrade without setting anything:**

- Read counts and depth drop, because by-product is removed. Nanopore keeps a median 55% of reads and
  72% of mean depth.
- The mapped span of a full-length product equals the product size in `CRM_Nested_products.bed`.
- Read names carry the amplicon their primer pair identified (`..._Amp7`).
- On Illumina, a pair is kept only if **both** mates pass filtering (no mate rescue). A library whose
  read 1 is a dark cluster over one amplicon can lose that amplicon; check the run summary.

**To keep the old behaviour:** `TRIM_MODE=standard`, or pick it in the launcher.

Validation: across 89 libraries, reported calls were identical to v1.2.1 in all 4 Illumina libraries
and in 79 of 81 real nanopore libraries (the 2 exceptions were low-yield flow-cell libraries that had
already failed, mean depth 72x and 20x). Amplicon sizes match the BED for all 10 amplicons; exact
spans are 98.7–98.9% (Illumina) and 96.3% (nanopore).

## 2. New and changed environment variables

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

## 3. Command-line changes, step by step

### Illumina (`MARK-I.sh`)

| step | v1.2.1 | v1.2.6 |
|---|---|---|
| merge | `fastp --merge` keeping merged reads, unmerged pairs and rescued single mates | same call, but in exact mode `--unpaired1/--unpaired2` are not used: both mates must pass |
| unmerged pairs | trimmed with `cutadapt -a/-A`, `-g/-G` and aligned as pairs | v1.2.3: not aligned at all. v1.2.4+ in exact mode: kept only when read 1 and read 2 start with the two primers of one amplicon (`cutadapt --pair-adapters --action=none`) |
| second merge | – | leftover pairs are 3'-quality-trimmed and merged again: `fastp --merge --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --disable_quality_filtering --length_required 30` |
| merged reads | `cutadapt -a file:<adapters>` then `-g file:<adapters>` | exact mode: `cutadapt -g file:<anchored primers> --revcomp --action=retain --rename '{id}_{adapter_name}'`, entries `^FORWARD...REVERSE$` |
| alignment | `bwa-mem2 mem` | same; in exact mode merged reads and primer-verified pairs are aligned separately and combined with `samtools merge` |

### Nanopore (`MARK.sh`)

| step | v1.2.1 | v1.2.6 |
|---|---|---|
| length safeguard | `fastplong --length_limit 1500` | same, plus `--disable_adapter_trimming` in exact mode (fastplong otherwise guesses an adapter and trims 10 bases past it, cutting into the primer; this also switches off its chimera splitting) |
| adapter / primer trimming | `cutadapt -b file:<adapters> --times 2` | exact mode: `cutadapt -g file:<linked primers> --revcomp --action=retain --rename '{id}_{adapter_name}' -e 0.15 -O 18`, then a per-read check that what lies outside each primer is library adapter and that no adapter sits inside the product |
| alignment | `minimap2 -ax map-ont` | same, plus `--end-bonus 15` in exact mode |

### Both

- v1.2.2: tool steps run in a plain shell (`bash -c`) instead of a login shell (`bash -lc`). A login
  shell re-read `/etc/profile`, which on Debian and Ubuntu reset `PATH` and hid the conda environment's
  `bwa-mem2`, `samtools` and `bcftools`.

## 4. New outputs

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

## 5. Launcher (`MARKLaunch.py`)

- `TRIM_MODE` is the first row of the parameter table, as an adapters/primers drop-down
  (`standard` / `exact`). Before this, the mode could not be chosen from the launcher.
- `FASTPLONG_ADAPTER_TRIM` is listed too (greyed out when the Illumina script is selected).
- Every parameter has a **"?"** beside its name. Clicking it explains what the setting does, which
  pipeline reads it, and its default in the selected script.

## 6. Tool versions

The conda recipe pins minimum versions to the set this release was validated with: fastp 1.3.7,
fastplong 0.7.1, cutadapt 5.2, bwa-mem2 2.3, minimap2 2.31, samtools 1.24, bcftools 1.24,
FastQC 0.12.1.

On macOS, cutadapt in multi-core mode can occasionally stop writing output part way through a large
nanopore file. The watchdog (`CUTADAPT_STALL_SECS`, default 600 s) catches this and reruns the step
single-core; the output is identical.

## 7. Known limitations

- **Depth is lower in exact mode**, by design: by-product is removed. Nanopore keeps a median 55% of
  reads and 72% of mean depth.
- **No mate rescue in exact mode.** If read 1 is a dark cluster over one amplicon, that amplicon can
  lose almost all of its depth even though read 2 is fine. Check the per-amplicon counts in the run
  summary, or use `TRIM_MODE=standard` for such a library.
- **fastplong's own quality filter** (its Q15 default, which MARK does not set) removes 3–4% of whole
  products on good nanopore libraries, and far more on low-quality ones.
- Reported calls are SNVs; indels are not reported in the main VCFs, unchanged from earlier releases.

## 8. Upgrading

- Nothing to change for calling: the variant filters are the same.
- Expect fewer reads and lower depth than v1.2.1 on the same data, because by-product is removed.
- To reproduce v1.2.1-style read retention, set `TRIM_MODE=standard`.
