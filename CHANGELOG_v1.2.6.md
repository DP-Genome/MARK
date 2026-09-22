# MARK v1.2.6 — fastplong adapter auto-trimming off in primers mode

This release includes everything in v1.2.5 (see `CHANGELOG_v1.2.5.md`). v1.2.5 remains available as a
fallback. Only nanopore runs with `TRIM_MODE=primers` change. Adapters mode and Illumina are
unchanged.

## What changed

By default fastplong guesses a library adapter for each file and trims it, plus 10 more bases. In
primers mode those 10 bases are the start of the primer, so the read then fails the primer check or
is kept clipped short. The primer step already removes everything outside the primers, so in primers
mode the auto-trimming adds nothing and costs good reads.

| | v1.2.5 | v1.2.6 |
|---|---|---|
| primers mode | fastplong adapter auto-trimming on | **off** (`--disable_adapter_trimming`) |
| adapters mode | on | on (unchanged) |
| new setting | – | `FASTPLONG_ADAPTER_TRIM` (default `0`); `1` turns the auto-trimming back on in primers mode |

- The setting is also in the launcher's parameter table, as `FASTPLONG_ADAPTER_TRIM`.
- The launcher can now choose the trimming mode. `TRIM_MODE` is the first row of its parameter table,
  as a drop-down (adapters / primers) that follows the script default (adapters) until overridden.
  Before this, runs started from the launcher were always in adapters mode.
- The run summary header records `TRIM_MODE` and whether fastplong auto-trimming was on.
- fastplong's chimera splitting is unchanged (on). Switching it off made no difference to real
  libraries in either mode.

## Evidence

The auto-trimming only acts where fastplong's guess finds an adapter: in the validation data, the 9
flongle HAC libraries and 3 unclassified bins. It touched 0 reads in every SUP library. With it off in
primers mode (flongle HAC 1–9):

- reads kept +82–100%, lowest amplicon about 2.5x
- reads not spanning one whole product 37–92 → 1–5 per library
- flongle_hac_1: 2,256 → 4,107 reads
- flongle_hac_3 recovers 16247A>G, 16261C>T, 16519T>C and 263A>G (the four calls it had lost against
  v1.2.1), plus 16217T>C
- error rate +0.1–0.3 points
- no other call changes

Details: `development/test_builds/MARK_FastplongDefaults_Test_20260922`.

## Verified

On flongle_hac_1 and MTC barcode16:

| case | flongle_hac_1 reads | MTC barcode16 reads |
|---|---|---|
| primers mode, default | 4,107 | 8,327 (unchanged) |
| primers mode, `FASTPLONG_ADAPTER_TRIM=1` | 2,256 (identical to v1.2.5) | 8,327 |
| adapters mode | 7,812 (identical to v1.2.5) | 13,042 |

## Also

- Version 1.2.6 in `MARK.sh`, `MARK-I.sh`, `MARKLaunch.py` and `conda-recipe/meta.yaml`.
