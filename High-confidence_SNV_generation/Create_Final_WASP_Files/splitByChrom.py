#!/usr/bin/env python

# Takes a snp file and formats it for WASP, splitting it into one file per scaffold/chromosome

# Usage: python3 splitScaffolds.py inFile

import os
import sys

# Read arguments
inFile = sys.argv[1]

# Create output file
outFileSuffix = ".snps.txt"

prevChrom = ""

outFile = "temp.txt"
outFile = open(outFile, 'w')


# Read input file
with open(inFile, 'r') as inFile:

	# Read first data line to initialize
	for line in inFile:
		chrom, p1, p2, refAlt = line.strip().split('\t')
		ref, alt = refAlt.split("|")
		
		outList = [p2, ref, alt]
		outPut = '\t'.join(outList)

		if chrom != prevChrom:
			outFile.close()
			outFile = chrom + outFileSuffix
			outFile = open(outFile, 'w')
			outFile.write(outPut + '\n')

		else: 
			outFile.write(outPut + '\n')

		prevChrom = chrom

outFile.close()

os.remove("temp.txt")
