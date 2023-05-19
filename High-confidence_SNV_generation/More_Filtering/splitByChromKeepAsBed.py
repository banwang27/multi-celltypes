# Usage: python3 splitScaffolds.py inFile species type_of_variant (indels or snps) type_of_comparison (e.g. huch, hugo, hurh etc. where huch indicates that the MAF was a human chimp alignment)

import os
import sys
import pandas as pd
import numpy as np

# Read arguments
inFile = sys.argv[1]
species = sys.argv[2]
suffix = sys.argv[3]

valid_chr = ["chr2A", "chr2B", "chrX", "chrY", "chrM"]
for i in range(23):
    valid_chr.append("chr" + str(i+1))
# Read input file
d = {}
with open(inFile, 'r') as inFile:

# Read first data line to initialize
    for line in inFile:
        chrom, p1, p2, refAlt = line.strip().split('\t')
        ref, alt = refAlt.split("|")
        if "chr" not in chrom:
            chrom = "chr" + str(chrom)
        outList = [chrom, int(p1), int(p2), ref, alt]
        if chrom in d.keys():
            d[chrom].append(outList)
        else:
            d[chrom] = [outList]
for key in d.keys():
    print(key)
    if key in valid_chr:
        df = pd.DataFrame(d[key])
        df.columns = ["chr", "start", "start2", "dont care", "bleh"]
        df = df.sort_values(by=['start'])
        outFile = open(key + "_" + species + suffix, "w+")
        for index, row in df.iterrows():
            line = row[0] + "\t" + str(row[1]) + "\t" + str(row[2]) + "\t" + row[3] + "|" + row[4]
            outFile.write(line)
            outFile.write("\n")
        outFile.close()