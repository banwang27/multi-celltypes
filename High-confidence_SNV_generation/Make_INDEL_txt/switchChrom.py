#!/usr/bin/env python

# Changes chromosome naming convention, chromosome must be in the first column

# Usage: python sbatch_all.py in out chroms mode

# Alllows us to run terminal commands and read arguments
import sys

inFile, outFile, chroms, mode = sys.argv[1:]

chromDict = {}

with open(chroms, 'r') as chroms: 
    for line in chroms:
        refseq, ens = line.strip().split('\t')
        if mode == 'refseq':
            chromDict[ens] = refseq
        elif mode == 'ensembl':
            chromDict[refseq] = ens


out = open(outFile, 'w')

with open(inFile, 'r') as inFile:
    for line in inFile:
        chrom = line.strip().split('\t')[0]
        fields = line.strip().split('\t')[1:]
        if chrom in chromDict.keys():
            newChrom = chromDict[chrom]

            outlist = [newChrom] + fields
            output = '\t'.join(outlist) + '\n'
            out.write(output)

out.close()
