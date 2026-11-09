#!/usr/bin/env python3
"""MARK v1.7 mixture check: is there a second haplotype in the library?

A library that is a mixture of two people (or carries contamination) has, besides its own sequence, reads from a
second mtDNA. The second person's bases at the positions where the two differ then sit TOGETHER on the same reads.
Sequencing errors and DNA damage also put minor bases at many positions, but each on its own read, so they are not
linked. This check looks for minor bases that are linked.

For the final (primer-clipped) BAM:
  1. At every position outside homopolymer regions and excluded sites, the most common base other than the majority
     base among bases of quality >= BASEQ_MIN (reads of mapping quality >= MAPQ_MIN) is a candidate when it is on at
     least 1% of those reads and on at least MIN_READS reads.
  2. For every two candidates that are covered together by at least MIN_SPAN reads: n_both = reads with the minor base
     at both. The pair is linked when n_both >= MIN_READS, n_both is at least LINK (70%) of the smaller of the two
     minor-base read counts on those reads (unlinked errors would share about AF1 x AF2 of them), and
     n_both / reads covering both >= MIX_MIN_AF (the share of reads carrying the second haplotype).
  3. The library is flagged when at least one pair is linked. When all linked pairs are in one amplicon the warning
     says so: that can be contamination concentrated there, or an artefact of that amplicon. The estimated share of
     the second haplotype is the median n_both / reads covering both over the linked pairs (in each amplicon it
     depends on that amplicon's depth in each sample, so it can differ between amplicons).

usage: MARK_mixture_check.py <final.bam> <regions.bed> <out.tsv> <MIX_MIN_AF> [BASEQ_MIN] [MAPQ_MIN]
Writes <out.tsv> (one line per linked pair, with a summary in its header lines) and prints one summary line:
  MIXTURE<TAB>yes|no<TAB>estimated share<TAB>n positions<TAB>positions<TAB>amplicons
"""
import collections, itertools, statistics, sys
import pysam

CAND_AF, MIN_READS, MIN_SPAN, LINK = 0.01, 5, 50, 0.7
real = lambda lin: (lin + 8284) % 16569 or 16569          # linearised reference -> rCRS


def main(bam_path, bed, out, mix_min_af, baseq=20, mapq=20):
    mask = set()
    for line in open(bed):
        c = line.split()
        if len(c) >= 4 and not line.startswith("#"):
            mask.update(range(int(c[1]) + 1, int(c[2]) + 1))
    bam = pysam.AlignmentFile(bam_path)
    count = collections.defaultdict(collections.Counter)
    reads = []                                             # per read: {position: base}, amplicon
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < mapq:
            continue
        q, s, d = a.query_qualities, a.query_sequence, {}
        for qp, rp in a.get_aligned_pairs(matches_only=True):
            p = rp + 1
            if p in mask or q[qp] < baseq:
                continue
            d[p] = s[qp]; count[p][s[qp]] += 1
        amp = a.query_name.rsplit("_", 1)[-1] if "_Amp" in a.query_name else ""
        reads.append((d, amp))
    cand = {}
    for p, c in count.items():
        n = sum(c.values()); top = c.most_common(2)
        if len(top) == 2 and top[1][1] >= MIN_READS and top[1][1] / n >= CAND_AF:
            cand[p] = top[1][0]
    pair = collections.defaultdict(lambda: [0, 0, 0, 0, collections.Counter()])   # span, minor1, minor2, both, amps
    order = sorted(cand)
    for d, amp in reads:
        hit = [p for p in order if p in d]
        for i, j in itertools.combinations(hit, 2):
            t = pair[(i, j)]; m1 = d[i] == cand[i]; m2 = d[j] == cand[j]
            t[0] += 1; t[1] += m1; t[2] += m2; t[3] += m1 and m2
            if m1 and m2: t[4][amp] += 1
    linked = []
    for (i, j), (span, m1, m2, both, amps) in pair.items():
        if span >= MIN_SPAN and both >= MIN_READS and both >= LINK * max(1, min(m1, m2)) and both / span >= mix_min_af:
            linked.append((i, j, span, m1, m2, both, amps.most_common(1)[0][0] if amps else ""))
    linked.sort()
    pos = sorted({x for i, j, *_ in linked for x in (i, j)})                 # panel order: HV1 first
    amps = sorted({r[6] for r in linked if r[6]}, key=lambda x: int(x[3:]) if x[3:].isdigit() else 0)
    share = statistics.median(b / s for i, j, s, m1, m2, b, a in linked) if linked else 0.0
    names = [f"{real(p)}{cand[p]}" for p in pos]
    with open(out, "w") as f:
        f.write(f"# MARK mixture check: {'EVIDENCE OF A SECOND HAPLOTYPE' if linked else 'no evidence of a second haplotype'}"
                f" (minimum share MIX_MIN_AF = {mix_min_af:g})\n")
        if linked:
            f.write(f"# estimated share of reads from the second haplotype: {share:.1%} (median over {len(linked)} linked pairs)\n"
                    f"# its bases: {' '.join(names)}\n# amplicons: {' '.join(amps) or '-'}\n")
            if len(amps) == 1:
                f.write(f"# seen in one amplicon only ({amps[0]}): a second haplotype in a whole library (a mixture or "
                        f"contamination) usually shows in several amplicons; in one amplicon only it may be contamination "
                        f"concentrated there, or an artefact of that amplicon (errors linked on its reads). Check {amps[0]}.\n")
        f.write("position_1\tbase_1\tposition_2\tbase_2\tamplicon\treads_covering_both\treads_minor_1\treads_minor_2"
                "\treads_minor_both\tshare\n")
        for i, j, span, m1, m2, both, amp in linked:
            f.write(f"{real(i)}\t{cand[i]}\t{real(j)}\t{cand[j]}\t{amp}\t{span}\t{m1}\t{m2}\t{both}\t{both / span:.4f}\n")
    print(f"MIXTURE\t{'yes' if linked else 'no'}\t{share:.4f}\t{len(pos)}\t{' '.join(names)}\t{' '.join(amps)}")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], a[2], float(a[3]), *(int(x) for x in a[4:6]))
