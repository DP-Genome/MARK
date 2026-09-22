# MARK v1.2.6 — clearer trimming-mode names; fastplong adapter auto-trimming off in exact mode

This release includes everything in v1.2.5 (see `CHANGELOG_v1.2.5.md`). v1.2.5 remains available as a
fallback. Only nanopore runs with `TRIM_MODE=exact` change. Standard mode and Illumina are
unchanged.

## Trimming modes renamed

Both modes remove the library adapters, and in both the primer footprints are soft-clipped before
calling, so "adapters" versus "primers" read like a choice between two things every run needs. The
modes differ in how a read is cut and which reads are kept, and are now named for that:

| v1.2.4–v1.2.5 | v1.2.6 | what it does |
|---|---|---|
| `TRIM_MODE=adapters` | `TRIM_MODE=standard` (default) | adapter sequences from the adapter list are cut wherever found; any read long enough afterwards is kept |
| `TRIM_MODE=primers` | `TRIM_MODE=exact` | exact amplicon trimming: a read is kept only with both primers of one amplicon, and cut exactly at them |

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
  as a drop-down (standard / exact) that follows the script default (standard) until overridden.
  Before this, runs started from the launcher were always in standard mode.
- The run summary header records `TRIM_MODE` and whether fastplong auto-trimming was on.
- fastplong's chimera splitting is unchanged (on). Switching it off made no difference to real
  libraries in either mode.

## Evidence

The auto-trimming only acts where fastplong's guess finds an adapter: in the validation data, the 9
flongle HAC libraries and 3 unclassified bins. It touched 0 reads in every SUP library. With it off in
exact mode (flongle HAC 1–9):

- reads kept +82–100%, lowest amplicon about 2.5x
- reads not spanning one whole product 37–92 → 1–5 per library
- flongle_hac_1: 2,256 → 4,107 reads
- flongle_hac_3 recovers 16247A>G, 16261C>T, 16519T>C and 263A>G (the four calls it had lost against
  v1.2.1), plus 16217T>C
- error rate +0.1–0.3 points
- no other call changes

Details: `development/test_builds/MARK_FastplongDefaults_Test_20260922`.

## Verified

- **Names:** `TRIM_MODE=exact` and the old `TRIM_MODE=primers` give identical results (MTC barcode16:
  8,327 reads). The default is standard (13,042). `TRIM_MODE=foo` stops with
  "TRIM_MODE must be standard or exact".
- **fastplong auto-trimming:**

On flongle_hac_1 and MTC barcode16:

| case | flongle_hac_1 reads | MTC barcode16 reads |
|---|---|---|
| exact mode, default | 4,107 | 8,327 (unchanged) |
| exact mode, `FASTPLONG_ADAPTER_TRIM=1` | 2,256 (identical to v1.2.5) | 8,327 |
| standard mode | 7,812 (identical to v1.2.5) | 13,042 |

## Also

- Version 1.2.6 in `MARK.sh`, `MARK-I.sh`, `MARKLaunch.py` and `conda-recipe/meta.yaml`.
