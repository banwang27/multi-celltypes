# Usage: python3 splitScaffolds.py inFile species type_of_variant (indels or snps) type_of_comparison (e.g. huch, hugo, hurh etc. where huch indicates that the MAF was a human chimp alignment)

import os
import sys
import pandas as pd

# Read arguments
inFile = sys.argv[1]
species = sys.argv[2]
type_of_variant = sys.argv[3]
type_of_comparison = sys.argv[4]

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
		outList = [p2, ref, alt]
		if chrom in d.keys():
		    d[chrom].append(outList)
		else:
		    d[chrom] = [outList]
for key in d.keys():
    print(key)
    if key in valid_chr:
        df = pd.DataFrame(d[key])
        df.to_csv(key + "_" + type_of_comparison + "_" + species + "." + type_of_variant +".txt", sep = "\t", index = False, header = False)
