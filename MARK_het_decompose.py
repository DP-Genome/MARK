#!/usr/bin/env python3
"""MARK v1.7 heteroplasmy report: split freebayes alleles into single events by their CIGAR.

freebayes reports a haplotype allele when events sit close together, e.g. at rCRS 66-73 the allele TGGGGGGTA>TGGGGGTG
(CIGAR 1M1D6M1X): one G fewer in the G6 tract and 73G. Splitting such an allele by comparing the strings base by base
(bcftools norm --atomize) shifts every base after the deletion and creates substitutions that no read carries (here
71T and 72G), and it leaves 73G in two pieces, so a homoplasmic 73G is reported as a 90% heteroplasmy. This script
splits each alternative allele by its CIGAR (INFO/CIGAR, one per allele): X gives substitutions, I an insertion and D
a deletion (both anchored on the base before them), M nothing. Every new record carries the reads of its allele
(FORMAT AD = reads with the reference, reads with the allele; GT 0/1, GL dropped, QA of that allele), so the reads of one event that sits in
several alleles are added up later. Alleles without a CIGAR are passed on unchanged.
For a substitution, FORMAT DP is set to the reads that have a base at that position: the reads of the reference and of
every allele that does not delete it (RO + those AO). Reads with a deletion there are a length variant (reported as
such), not a second base, so they are not counted against it; e.g. DNA 007 460C sits in a C8 tract, and reads with one
C fewer would otherwise turn a homoplasmic 460C into a 70% heteroplasmy.
usage: freebayes ... | MARK_het_decompose.py > out.vcf     (then bcftools norm -m -any -f ref to left-align indels)
"""
import sys

def deleted(cigar, pos):
    """reference positions an allele deletes"""
    out, p = set(), pos
    for n, op in ((int(c[:-1]), c[-1]) for c in __import__("re").findall(r"\d+[MXID]", cigar)):
        if op in "MX": p += n
        elif op == "D": out.update(range(p, p + n)); p += n
    return out

def split(ref, alt, cigar, pos):
    out, r, a, p = [], 0, 0, pos
    for n, op in ((int(c[:-1]), c[-1]) for c in __import__("re").findall(r"\d+[MXID]", cigar)):
        if op == "M":
            r += n; a += n; p += n
        elif op == "X":
            for k in range(n):
                out.append((p + k, ref[r + k], alt[a + k]))
            r += n; a += n; p += n
        elif op == "I":
            anc = ref[r - 1] if r else ref[0]
            out.append((p - 1, anc, anc + alt[a:a + n]) if r else (p, ref[0], alt[a:a + n] + ref[0])); a += n
        elif op == "D":
            anc = ref[r - 1] if r else None
            out.append((p - 1, anc + ref[r:r + n], anc) if r else (p, ref[r:r + n + 1], ref[r + n])); r += n; p += n
    return out

for line in sys.stdin:
    if line.startswith("#"):
        sys.stdout.write(line); continue
    x = line.rstrip("\n").split("\t")
    info = dict(kv.split("=", 1) for kv in x[7].split(";") if "=" in kv)
    alts = x[4].split(","); cig = info.get("CIGAR", "").split(",")
    fk, fv = x[8].split(":"), x[9].split(":")
    fmt = dict(zip(fk, fv)); ad = fmt.get("AD", ".").split(",")
    if len(cig) != len(alts) or len(ad) != len(alts) + 1:
        sys.stdout.write(line); continue
    dels = [deleted(c, int(x[1])) for c in cig]
    for i, alt in enumerate(alts):
        for p, rr, aa in split(x[3], alt, cig[i], int(x[1])):
            y = list(x); y[1], y[3], y[4] = str(p), rr, aa
            g = dict(fmt); g["AD"] = f"{ad[0]},{ad[i + 1]}"
            if len(rr) == len(aa) == 1 and "DP" in g and all(v.isdigit() for v in ad):
                g["DP"] = str(int(ad[0]) + sum(int(ad[j + 1]) for j in range(len(alts)) if p not in dels[j]))
            if "AO" in g: g["AO"] = ad[i + 1]
            if "GT" in g: g["GT"] = "0/1"                      # one alternative allele per record now
            if "GL" in g: g["GL"] = "."                        # likelihoods of the multi-allele record do not apply
            for k in ("QA",):
                if k in g and "," in g[k]: g[k] = g[k].split(",")[i]
            y[8] = ":".join(fk); y[9] = ":".join(g[k] for k in fk)
            sys.stdout.write("\t".join(y) + "\n")
