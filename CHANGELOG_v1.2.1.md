# v1.2.1 — primer-match tolerance raised from 6 to 10

## What changed

One number, in both pipelines:

```
PRIMER_TOL = 6   ->   PRIMER_TOL = 10
```

Nothing else. The insert table, the `products` table, both primer boundaries per amplicon, the
adapter lists and every processing step are unchanged from v1.2.

## Why

`PRIMER_TOL` is how close a read end must be to a known PCR product terminus before the trimmer
treats that end as carrying the product's primer. At 6, one population of reads was missed.

Amplicon 6's forward primer begins at lin **8249** in the untrimmed molecule. Cutadapt removes
its first 9 bases, so **after adapter trimming — which is what the read trimmer actually sees —
the read starts at 8258**, and 8258 is what the `products` table correctly records. Measured in
2800M's pre-trim alignment: 85,474 reads start at 8258.

But **1,525 reads escape adapter trimming** and still start at 8249. Those are 9 bases away from
the table entry, past a tolerance of 6, so they matched nothing, fell through to the max-overlap
fallback, and kept amplicon 6's primer bases at lin 8249-8270 (rCRS 16,533-16,553). In the final
v1.2 BAM, 1,197 such reads survive with primer attached.

Amplicon 9 has the same structure — raw start 8602, adapter-trimmed start 8608, a 6-base shift —
and matched only because 6 was exactly equal to the tolerance, with zero margin.

Raising the tolerance to 10 catches both populations of both amplicons.

## The size of the effect depends on the platform

The paragraphs above describe Illumina, where only 1.7% of amplicon 6 reads escape adapter
trimming. **On nanopore it is the majority.** The nanopore adapter trimmer does not remove the
first nine bases of amplicon 6's primer, so 81-87% of reads start at the true primer position,
8249, rather than 8258. At a tolerance of 6 they all missed the table entry and kept their primer.

Measured across all 85 libraries, the share of v1.2's depth at rCRS 16,534-16,554 that was
primer: **73% on MinION / MTC libraries (54), 77% on flongle (27), 3% on Illumina (4).** In
v1.2.1 the depth across that footprint matches the flanking sequence within 0.4%.

Consequences across all 85 libraries: mean depth −4.9%, amplicon 5 −17% (the footprint sits
inside amplicon 5's insert), imbalance 5.3x → 4.4x, variant calls identical in 85/85. The
earlier "−0.6%" figure came from 13 libraries of which only nine were nanopore, and understated
the effect.

## Why 10, and why it is safe

The closest two product starts are 79 bases apart, the closest two product ends 71 bases apart,
and the closest start-to-end pair 13 bases apart. A tolerance of 10 therefore cannot make one
amplicon's terminus match another's, nor a start match an end. It has margin at both ends: the
largest real shift is 9 (amplicon 6), and the smallest distance to a wrong match is 13.

## A correction to an earlier conclusion

An earlier draft of this release changed the `products` table instead, moving Amp6 to 8249 and
Amp9 to 8602, on the basis that those are the starts observed in a raw alignment. **That was
wrong and has been reverted.** The raw alignment shows the molecule; the trimmer operates on
adapter-trimmed reads, and the table has to describe what the trimmer sees. Making that change
broke the match for the 85,474 reads that legitimately start at 8258 and left 5,359 of them
carrying primer — measured as roughly 2.6x more retained primer sequence than v1.2, not less.

The lesson worth keeping: **coordinates in the `products` table are post-adapter-trim
coordinates.** Any future re-derivation must be done on the pipeline's own `_initial_sorted.bam`,
never on a raw alignment.

## Verification before reprocessing

| fragment | in | TOL=6 | TOL=10 | correct |
|---|---|---|---|---|
| Amp6, adapter-trimmed start | 8258-8465 | 8271-8437 | 8271-8437 | 8271-8437 |
| **Amp6, adapter escaped** | 8249-8465 | **8249-8437** | **8271-8437** | 8271-8437 |
| **Amp6-F escaped + Amp5-R** | 8249-8337 | **8249-8315** | **8271-8315** | 8271-8315 |
| Amp6-F trimmed + Amp5-R | 8258-8337 | 8271-8315 | 8271-8315 | 8271-8315 |
| Amp9, adapter-trimmed start | 8608-8745 | 8627-8721 | 8627-8721 | 8627-8721 |
| Amp9, adapter escaped | 8602-8745 | 8627-8721 | 8627-8721 | 8627-8721 |
| Amp4-F + Amp5-R hybrid | 8079-8337 | 8103-8315 | 8103-8315 | 8103-8315 |
| Amp10 truncated forward read | 8687-8723 | 8714-8723 | 8714-8723 | 8714-8723 |
| read inside Amp9's reverse primer | 8727-8745 | dropped | dropped | dropped |

The two rows in bold are what this release fixes.

## What the published coordinates confirm

Vinueza-Espinosa et al. 2023 (*Electrophoresis* 44:1423-1434) publishes this kit's ten amplicons
in rCRS. Converted with `lin = (rCRS - 8284) mod 16569`, **all ten match our insert table
exactly**, including the two that wrap the origin (Amp5 16474-30, Amp6 16555-152). That is an
independent confirmation that the published coordinates are the inserts, which this pipeline has
always assumed.

The paper publishes inserts only, and Promega does not release the CRM primer sequences, so
product spans still have to be derived from the data. The paper's stated amplicon size range of
147-237 bp is a useful cross-check: eight of ours fall inside it, with Amp3 at the 237 maximum
and Amp4 at the 147 minimum.

On amplicon performance the paper reports amplicons **2, 3 and 8** as weakest — 3 because it is
the largest and worst affected by fragmentation, 2 and 8 because they sit on the damage and
heteroplasmy hotspots at 16189 and 303-315. **Amplicon 10 is not flagged**, supporting our own
finding that DNA007's amplicon 10 failure is specific to that sample rather than a kit weakness.

## Also in v1.2.1 — cutadapt can no longer hang a run

cutadapt's multi-core mode can deadlock on macOS: its parent process waits forever on worker
processes that have already died. It happened once during v1.2 testing — one library sat for
four days having used four seconds of CPU — and because the pipeline had no timeout, the run
simply stopped, with no error.

Every cutadapt call now goes through `run_cutadapt`, which watches the output file. If the
output stops growing for `CUTADAPT_STALL_SECS` (default 600 s), the run is killed and repeated
with `--cores 1`, which has no worker processes and cannot deadlock. If it stalls again, the
pipeline stops with an error rather than skipping the sample.

The retry cannot change results, because cutadapt keeps reads in input order in multi-core
mode. On `Test_M.fastq`, output through the wrapper, direct at 8 cores and direct at 1 core
are byte-identical (4,846 reads, same MD5). Also tested: a simulated deadlock (killed, rerun
single-core, run completes), a deadlock that persists single-core (pipeline stops with an
error), and a cutadapt failure (exit code passed through, pipeline stops).

## Also in v1.2.1 — the scripts use their own data files, from any directory

Up to v1.2 the reference, regions BED and adapter list were looked up in the current working
directory, or beside the input folder. Started from anywhere else, a command-line run either
stopped — `Error: 'linearized_regions.bed' not found` — or silently used whatever copy happened
to be in that folder, which is exactly how a stale adapter list could creep back in.

Bare file names now resolve to the copy shipped beside the script first, following symlinks,
so a conda install or a repository checkout uses its own files from any directory. Explicit
settings still win: `ref=`, `regions_bed=` and `ADAPTER_FILE=` given as paths are used exactly
as given. The dashboard, which already passes all three explicitly, is unaffected, and its
parser still reads the same defaults from the scripts. A custom adapter list must now be given
through `ADAPTER_FILE` (or the dashboard's adapter field) rather than by placing a file named
`MARK_Adapter_List_*.txt` next to the input.

The resolved path of each file is recorded in the run summary under `## Reference`,
`## Regions BED` and `## Adapter File`.

Tested from an empty folder with no settings: before the fix the installed `MARK.sh` exited with
the error above; after it, `MARK.sh` and `MARK-I.sh` both complete, resolve all three files beside
the script, write nothing into the working directory, and still honour an explicit
`ADAPTER_FILE`.
