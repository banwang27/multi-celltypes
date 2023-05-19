# usr/bin/env python

# Filter SNP lists for SNps that we want to keep

# usage python SNPfilter.py infile keep_list outfile

import sys

infile, keep, outfile = sys.argv[1:]

snp_dict = {}
print("Reading SNPs...")
with open(infile, 'r') as infile:
    for line in infile:
        chrom, pos0, pos1, snp = line.strip().split()
        if chrom in snp_dict:
            snp_dict[chrom][pos1] = snp
        else:
            snp_dict[chrom] = {pos1: snp}

print("Finished.")
print("Filtering...")

out = open(outfile, 'w')

with open(keep, 'r') as keepfile:
    for line in keepfile:
        chrom, pos = line.strip().split()
        if chrom in snp_dict:
            if pos in snp_dict[chrom]:
                output = "\t".join([chrom, str(int(pos)-1), pos, snp_dict[chrom][pos]])
                out.write(output + '\n')

print("Done.")
out.close()
