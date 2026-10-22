# MARK v1.4.4

**Launcher layout only: the Activity Log is a compact strip, so the whole parameter table fits on screen. The pipeline
and its results are identical to v1.4.3.**

## What changed
- **Activity Log:** a fixed strip about 170 px tall (it scrolls) at the foot of the dashboard, instead of a panel that
  stretched with the window and took about a quarter of it (290 px on a maximized MacBook Pro window).
- **Parameter table:** the tabs take all the spare height. With the window maximized, every parameter row, down to
  FREEBAYES, is visible without scrolling.
- **README:** the dashboard screenshot shows the current launcher (it was from v1.1.5).
- **Version:** 1.4.4 in MARK.sh, MARK-I.sh, the launcher title and header, and the conda recipe.

## Unchanged
MARK.sh and MARK-I.sh differ from v1.4.3 only in the version string. Trimming, clipping, variant calling, filters,
reports and output files are the same.

## Validation
- **Layout:** on a maximized 2056 × 1188 pt window, the Activity Log is 168 pt and the parameter table needs 823 of the
  836 pt available, so it does not scroll (v1.4.3: log 290 pt, table area 765 pt, table scrolled).
- **Pipeline:** `bash -n` passes for MARK.sh and MARK-I.sh; their diff against v1.4.3 is the version line.
