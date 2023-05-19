import pysam
import argparse
from tqdm import tqdm
from collections import defaultdict


def pipesplit(col):
    return lambda input: input.split("|")[col]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("snps", help="SNP BED file")
    parser.add_argument("reads", help="Mapped reads file BAM or (untested) SAM")
    parser.add_argument("output")
    return parser.parse_args()


def parse_bed(fname):
    outdict = defaultdict(dict)
    for line in open(fname):
        chr, pos0, pos1, refalt = line.split()
        pos0 = int(pos0)
        chrdict = outdict[chr]
        chrdict[pos0] = refalt.split("|")
    return outdict


if __name__ == "__main__":
    args = parse_args()
    snps = parse_bed(args.snps)

    out_table = {(chr, pos + 1): [0, 0, 0] for chr in snps for pos in snps[chr]}

    reads = pysam.AlignmentFile(args.reads)
    chrnames = reads.references
    for read in tqdm(reads, total=reads.mapped):
        chrom = chrnames[read.reference_id]
        if chrom not in snps:
            continue
        chrsnps = snps[chrom]
        seq = read.seq

        for rpos, pos, base in read.get_aligned_pairs(matches_only=True, with_seq=True):
            base = read.seq[rpos].upper()
            ptuple = chrom, pos + 1
            if pos in chrsnps:
                ref, alt = chrsnps[pos]
                if base == ref:
                    out_table[ptuple][0] += 1
                elif base == alt:
                    out_table[ptuple][1] += 1
                else:
                    out_table[ptuple][2] += 1

    with open(args.output, "w") as outfh:
        print("CHROM", "POS", "REF", "ALT", "NON_REFALT", sep="\t", file=outfh)
        for (chrom, pos) in sorted(out_table):
            print(chrom, pos, *out_table[(chrom, pos)], sep="\t", file=outfh)
