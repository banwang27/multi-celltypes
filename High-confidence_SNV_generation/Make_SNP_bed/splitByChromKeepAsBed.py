# Usage: python3 splitScaffolds.py inFile species type_of_variant (indels or snps) type_of_comparison (e.g. huch, hugo, hurh etc. where huch indicates that the MAF was a human chimp alignment)

import os
import sys
import pandas as pd

# Read arguments
inFile = sys.argv[1]
species = sys.argv[2]

valid_chr = ["chr2A", "chr2B", "chrX", "chrY", "chrM"]
for i in range(23):
    valid_chr.append("chr" + str(i+1))
# Read input file
d = {}
with open(inFile, 'r') as inFile:

	# Read first data line to initialize
	for line in inFile:
		chrom, p1, p2, refAlt = line.strip().split('\t')
		if "chr" not in chrom:
		    chrom = "chr" + str(chrom)
		outList = [chrom, p1, p2, ref, alt]
		if chrom in d.keys():
		    d[chrom].append(outList)
		else:
		    d[chrom] = [outList]
for key in d.keys():
    print(key)
    if key in valid_chr:
        df = pd.DataFrame(d[key])
        df = df.sort_values(by = 1)
        outFile = open("", "w+")
        for index, row in df.iterrows():
        	line = row[0] + "\t" + row[1] + "\t" + row[2] + "\t" + row[3] + "|" + row[4]
        	outFile.write(line)
            outFile.write("\n")
        outFile.close()