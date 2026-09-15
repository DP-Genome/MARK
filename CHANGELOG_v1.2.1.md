# MARK v1.2.1 — changes since v1.1.8

This release also includes the changes made in v1.1.9. The variant-calling settings and filters are unchanged. What changed is how reads are filtered, trimmed and aligned before calling: more real template reaches the caller, and adapter and primer sequence is no longer counted as sample.

## Adapter lists (both pipelines)

| | v1.1.8 | v1.2.1 |
|---|---|---|
| Illumina (`MARK-I.sh`) | `Updated_Adapter_Primer_List_Cutadapt_cleaned.txt`, 744 entries, shared | `MARK_Adapter_List_Illumina.txt`, 200 entries |
| Nanopore (`MARK.sh`) | the same shared list | `MARK_Adapter_List_ONT.txt`, 207 entries (Illumina + Nanopore adapters) |

- **Mitochondrial sequence removed.** The old list held 320 entries of mtDNA control-region sequence, plus several adapter–mtDNA chimeras. cutadapt `-a` removes a match *and everything after it*, and `-g` removes a match *and everything before it*. So each mtDNA entry cut every read containing it at a fixed position in the genome, and the remainder often failed the length filter. No entry now has a ≥15 bp exact match to the reference. `MARK_Adapter_List_v2.rebuild_log.txt` records every edit.
- **Index barcodes removed.** These 8 bp sample indices are read in a separate index cycle and never occur in the insert. With `--overlap 5`, a chance 5 bp match could trim a real read end.
- **Nanopore chemistry removed from the Illumina list.** Illumina libraries cannot contain it.
- The dashboard now switches to the matching list when you change pipelines, unless you have chosen your own.

## Nanopore adapter trimming (`MARK.sh`)

```
v1.1.8   cutadapt -a file:LIST ...   then   cutadapt -g file:LIST ...
v1.2.1   cutadapt -b file:LIST --times 2 ...
```

These amplicons carry adapter at both ends of the read. `-a` treats every match as a 3′ adapter, so an adapter matched at the start of a read removed the entire read, which `--minimum-length` then discarded. `-b` treats a match that includes the read's first base as a 5′ adapter, keeping the sequence after it, and any other match as a 3′ adapter. `--times 2` lets one read lose an adapter from each end. Length floors are unchanged (`MIN_LEN=90`, `MIN_LEN_POST=90`).

## Illumina read handling (`MARK-I.sh`)

```
v1.1.8   fastp --merge --include_unmerged --merged_out merged.fq --length_required 90 ...
v1.2.1   fastp --merge --trim_poly_g --merged_out merged.fq \
               --out1 unmerged_R1.fq --out2 unmerged_R2.fq \
               --unpaired1 rescued_R1.fq --unpaired2 rescued_R2.fq --length_required 90 ...
```

- **`--trim_poly_g`** removes poly-G tails, which two-colour Illumina chemistry produces when a cluster gives no signal. Left in, they align with most of the read soft-clipped and add depth that is not sample.
- **Unmerged pairs stay paired.** `--include_unmerged` placed them in the merged file as two unrelated single-end reads, so where mates overlapped, one fragment was counted twice. Now each unmerged pair goes through four steps:
  1. it is written to `--out1`/`--out2`;
  2. it is trimmed in paired mode (`cutadapt -a/-A`, then `-g/-G`);
  3. it is aligned as a pair (`bwa-mem2 mem ref R1 R2`);
  4. it is combined with the merged reads (`samtools merge`).

  After primer trimming, the BAM passes through `samtools sort -n | samtools fixmate | samtools sort`, so mate positions stay correct and the pileup counts each fragment once.
- **A good mate is no longer discarded with a bad one.** By default fastp drops the whole pair when either read fails filtering. `--unpaired1`/`--unpaired2` keep the read that passed, and it continues as single-end. This matters when a fault affects only one read of the pair, such as poly-G in R1.
- **The 90 bp floor applies once, to the raw fragment.** v1.1.8 reapplied it at every cutadapt step and at the post-trim filter, which discarded valid fragments once their adapter and primer were removed. After trimming, the floor is now `MIN_LEN_POST=30`, a minimum for reliable mapping (default was 90).

## Primer trimming (both pipelines)

**v1.1.8:** each read was assigned to the amplicon it overlapped most, then cut to fixed, non-overlapping tiles.

**v1.2.1:** each read end is compared with the known PCR product ends, within `PRIMER_TOL = 10` bp.
- An end at a product start has that product's forward primer removed.
- An end at a product end has its reverse primer removed.
- Everything between the primers is kept, whichever amplicons it spans.
- A read with neither end at a known product end falls back to the old rule and is clipped to a single insert.

This gives four improvements:

- **Fragments spanning two amplicons keep all their template.** Examples are a forward primer paired with the neighbouring amplicon's reverse primer. Previously, such a fragment was clipped to one amplicon.
- **Primer is removed by position, not by amplicon assignment.** Primer sequence always matches the reference, so counting it as coverage can hide a real variant under its footprint.
- **The insert table matches the ten inserts published for the PowerSeq CRM Nested System.** Adjacent inserts therefore overlap as they do in the kit, and a base covered by two amplicons counts reads from both.
- **Reads outside every amplicon are dropped as off-target.** Previously they crashed the trimming step.

The tolerance matters because the number of primer bases left at the start of a read after adapter trimming varies, and differs between Illumina and nanopore. The closest termini of two different amplicons are 13 bp apart, so a tolerance of 10 cannot match the wrong amplicon.

### Amplicon coordinates

The inserts are those published by Vinueza-Espinosa et al. 2023 (*Electrophoresis* 44:1423-1434). Promega does not publish the primer sequences, so the product spans are derived from the fixed read termini seen when sequencing this kit. The pipelines run against a linearized reference, `NC_012920.1_linearized`, where `rCRS = (lin + 8284) mod 16569`. Coordinates below are 1-based and inclusive.

| Amplicon | Insert (rCRS) | Insert (linearized) | PCR product (linearized) |
|---|---|---|---|
| Amp1 | 16013–16126 | 7729–7842 | 7702–7868 |
| Amp2 | 16116–16225 | 7832–7941 | 7810–7964 |
| Amp3 | 16223–16408 | 7939–8124 | 7913–8149 |
| Amp4 | 16387–16486 | 8103–8202 | 8079–8225 |
| Amp5 | 16474–30 | 8190–8315 | 8166–8337 |
| Amp6 | 16555–152 | 8271–8437 | 8249–8465 |
| Amp7 | 136–257 | 8421–8542 | 8394–8578 |
| Amp8 | 246–364 | 8531–8649 | 8503–8674 |
| Amp9 | 342–436 | 8627–8721 | 8602–8745 |
| Amp10 | 429–592 | 8714–8877 | 8687–8904 |

The same coordinates ship as BED files (0-based, half-open): `CRM_Nested_inserts.bed`, `CRM_Nested_products.bed` and `CRM_Nested_primers_empirical.bed`. The pipelines carry their own copy internally and do not read these files.

For Amp6 and Amp9, adapter trimming usually removes the first few primer bases, so reads start at 8258 and 8608 rather than at the product starts above. The trimmer's internal table records those positions, and `PRIMER_TOL` covers both cases.

### How the trimmer handles each kind of read

Linearized coordinates, 1-based.

| Read | What it is | Kept |
|---|---|---|
| 8258–8465 | Amp6 product, primer partly removed by adapter trimming | 8271–8437 |
| 8249–8465 | Amp6 product, full primer still attached | 8271–8437 |
| 8079–8337 | Amp4 forward primer to Amp5 reverse primer (hybrid) | 8103–8315 |
| 8249–8337 | Amp6 forward primer to Amp5 reverse primer (hybrid) | 8271–8315 |
| 8687–8723 | Amp10 read that stops early | 8714–8723 |
| 8727–8745 | read lying entirely inside Amp9's reverse primer | dropped |

## Variant calling (both pipelines)

`bcftools mpileup` now also receives `--max-idepth $PILEUP_MAX_IDEPTH`, which defaults to `PILEUP_MAX_DEPTH`. bcftools applies a separate depth cap of 250 to indel candidates. Raising `-d` alone left that cap in place, so no indels were proposed above 250×. This affects only the `annotated_all` and `qual_filtered` review VCFs; the `snps` and `clean` outputs remain SNP-only.

## Running the pipeline

- **cutadapt can no longer hang a run.** Its multi-core mode can deadlock on macOS. Every cutadapt call now watches its output file. If the output stops growing for `CUTADAPT_STALL_SECS` seconds, the call is rerun with `--cores 1`, and a second stall stops the pipeline with an error. Output is identical either way, because cutadapt keeps reads in input order.
- **The scripts use their own data files.** The reference, regions BED and adapter list are found beside the script, following symlinks, so command-line runs work from any directory. Explicit `ref=`, `regions_bed=` and `ADAPTER_FILE=` values still take priority. The resolved paths are recorded in `run_summary.txt`.
- **The run folder is created beside the input.** Previously a command-line run put it in the current directory, and a folder browsed in the dashboard could place it inside the input. `OUTPUT_DIR` sets another location, which is created if missing. `RUN_NAME` names the folder, and the dashboard's Custom Run Name field now works. A run will not write into an existing folder that is not empty.

## New and changed settings

| Setting | Default | Pipeline |
|---|---|---|
| `CUTADAPT_STALL_SECS` | 600 | both — new |
| `PILEUP_MAX_IDEPTH` | same as `PILEUP_MAX_DEPTH` | both — new |
| `OUTPUT_DIR` | the folder holding the input | both — new |
| `RUN_NAME` | auto-generated | both — new |
| `MIN_LEN_POST` | 30 (was 90) | `MARK-I.sh` |

## Upgrading from v1.1.8

- **Custom adapter lists:** pass them with `ADAPTER_FILE=` or the dashboard's adapter field. A list placed next to the input is no longer picked up, and the old file name is no longer searched for.
- **Output location:** run folders now appear beside the input, not in the directory the script was started from.
- **`run_summary.txt`:**
  - Illumina reports `2_Trim3(merged)` and `2b_TrimPairs` separately, plus a count of rescued mates.
  - Nanopore reports a single `2_AdapterTrim` step.
  - Both record `PILEUP_IDEPTH` and the resolved reference, regions BED and adapter list paths.

Version 1.2.1 across the dashboard, both pipelines and the conda recipe.
