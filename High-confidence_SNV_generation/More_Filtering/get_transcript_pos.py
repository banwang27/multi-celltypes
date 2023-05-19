import os
import sys
import pandas as pd

# Read arguments
inFile = sys.argv[1]
d = {}
species = inFile.split(".")[0]
chrs = ["chr2B", "chr2A", "chrX", "chrY"]
for i in range(1,23):
    chrs.append("chr" + str(i))
with open(inFile, 'r') as inFile:
    # Read first data line to initialize
    
    for line in inFile:
        if "#" != line[0]:
            chrom = line.strip().split('\t')[0]
            feat = line.strip().split('\t')[2]
            start = line.strip().split('\t')[3]
            end = line.strip().split('\t')[4]
            if feat == "transcript":
                if chrom in d.keys():
                    d[chrom].append([start, end])
                else:
                    d[chrom] = [[start, end]]
for key in d.keys():
    if key in chrs:
        outFile = key + "_" + species + "_transcript_positions.txt"
        df = pd.DataFrame(d[key]).drop_duplicates()
        df.to_csv(outFile, sep = "\t", index = False, header = False)