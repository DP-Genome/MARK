# MARK v1.2.6 — clearer trimming-mode names; fastplong adapter auto-trimming off in exact mode

This release includes everything in v1.2.5 (see `CHANGELOG_v1.2.5.md`). Two things change what a run
does by default: exact amplicon trimming is now the default mode, and fastplong's own adapter
auto-trimming is switched off in that mode. Variant-calling settings and filters are unchanged, and
`TRIM_MODE=standard` restores the pre-v1.2.4 behaviour.

## Exact amplicon trimming is now the default

`TRIM_MODE` defaults to `exact` from this release (it was `adapters`/`standard` in v1.2.4 and v1.2.5).
Exact trimming was validated against the v1.2.1 release on 89 libraries with identical calls, and it
is what removes PCR by-product from the BAM. To keep the pre-v1.2.4 behaviour, set
`TRIM_MODE=standard`, or choose it in the launcher.

What changes for a run that does not set `TRIM_MODE`:

- Reads are cut exactly at their own amplicon primers, and a read is kept only if it carries both
  primers of one amplicon. PCR by-product (primer dimers, read-through and amplicon-overlap products)
  no longer reaches the BAM, so **read counts and depth drop**: nanopore keeps a median 55% of reads
  and 72% of mean depth.
- Illumina pairs are kept only when both mates pass filtering (no mate rescue). A library whose read 1
  is a dark cluster over one amplicon can therefore lose that amplicon; check the run summary.
- Read names carry the amplicon their primer pair identified (`..._Amp7`).

## Trimming modes renamed

Both modes remove the library adapters, and in both the primer footprints are soft-clipped before
calling, so "adapters" versus "primers" read like a choice between two things every run needs. The
modes differ in how a read is cut and which reads are kept, and are now named for that:

| v1.2.4–v1.2.5 | v1.2.6 | what it does |
|---|---|---|
| `TRIM_MODE=adapters` | `TRIM_MODE=standard` | adapter sequences from the adapter list are cut wherever found; any read long enough afterwards is kept |
| `TRIM_MODE=primers` | `TRIM_MODE=exact` (**now the default**) | exact amplicon trimming: a read is kept only with both primers of one amplicon, and cut exactly at them |

The old names still work and print a one-line note. Results are unchanged. The run summary header
now records the mode with a short description, for Illumina as well as nanopore.

## fastplong adapter auto-trimming off in exact mode

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
  overshoot. It does not help: on the reads both settings keep, exact spans fall from 91.7% to 80.7%
  on low-quality libraries, and the ~250 reads it adds per library are only 55-58% exact while the
  ~180 it drops were 87-92% exact. fastplong's adapter detection is not precise enough on noisy
  reads to stop cleanly at the adapter. Details: development/test_builds/MARK_TrimExt_Test_20260923.

## Launcher: help for every parameter

Each parameter in the launcher's table now has a "?" beside its name. Clicking it opens a short card
that says what the setting does, which pipeline uses it (nanopore, Illumina or both), and its default
in the selected script. The card closes with Close, Escape, or a click elsewhere.
`DISCARD_WARN_PCT` is defined in both scripts but read by no step; its card says so.

## Evidence

The auto-trimming only acts where fastplong's guess finds an adapter: in the validation data, the 9
libraries basecalled with the fast model, and in unclassified barcode bins. It touched 0 reads in every
library basecalled with the high-accuracy model. With it off in exact mode, on 9 such libraries:

- reads kept +82–100%, lowest amplicon about 2.5x
- reads not spanning one whole product 37–92 → 1–5 per library
- one library went from 2,256 to 4,107 reads
- one library recovers the four calls it had lost against v1.2.1, plus one more
- error rate +0.1–0.3 points
- no other call changes

Details: `development/test_builds/MARK_FastplongDefaults_Test_20260922`.

## Verified

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

## Also

- Version 1.2.6 in `MARK.sh`, `MARK-I.sh`, `MARKLaunch.py` and `conda-recipe/meta.yaml`.
