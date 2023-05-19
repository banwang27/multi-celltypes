# Usage: python3 splitScaffolds.py inFile species type_of_variant (indels or snps) type_of_comparison (e.g. huch, hugo, hurh etc. where huch indicates that the MAF was a human chimp alignment)

import os
import sys
import pandas as pd

# Read arguments
Suffix1 = sys.argv[1]
Suffix2 = sys.argv[2]
species = sys.argv[3]
t = sys.argv[4]
min_reads = sys.argv[5]
if t == "transcript":
    outFile = open("GRCH38.SNPs.TRANSCRIPT.C3649.WGS.bed", "w+")
elif t == "exon":
    outFile = open("GRCH38.SNPs.EXON.C3649.WGS.bed", "w+")
elif t == "reads":
    outFile = open("GRCH38.SNPs." + "MIN_" + str(min_reads) + "_READS.C3649.WGS.bed", "w+")

chrs = []
if species == "chimp":
    chrs = ["chr2A", "chr2B", "chrX", "chrY"]
    for i in range(1, 23):
        if i != 2:
            chrs.append("chr" + str(i))
elif species == "human":
    chrs = ["chrX", "chrY"]
    for i in range(1, 23):
        chrs.append("chr" + str(i))

for chrom in chrs:
    bed = open(chrom + "_" + species + Suffix2 + ".bed", "r")
    l = []
    if t == "exon" or t == "transcript":
        text = open(chrom + "_" + species + "_" + Suffix1 + ".txt", "r")
        for line in text:
            try:
                start, end = line.split("\t")
                l.append([start, end])
            except:
                break
        cur_pos = 0
        for line in bed:
            if len(line.split("\t")) >= 2 and len(l[cur_pos]) >= 2:
                start = int(l[cur_pos][0])
                end = int(l[cur_pos][1])
                n = int(line.split("\t")[1])
                n_1 = int(line.split("\t")[2])
                while n_1 > end:
                    cur_pos += 1
                    start = int(l[cur_pos][0])
                    end = int(l[cur_pos][1])
                if n_1 >= start and n <= end:
                    outFile.write(line)
                else:
                    pass
    else:
        text = open(chrom + "_" + species + "MIN_" + min_reads + "_READS" + "_" + "KEEP.txt", "r")
        for line in text:
            l.append(line)
        cur_pos = 0
        for line in bed:
            if len(line.split("\t")) >= 2:
                if cur_pos >= len(l):
                        break
                else:
                    snp_pos = int(l[cur_pos])
                n = int(line.split("\t")[1])
                n_1 = int(line.split("\t")[2])
                while n_1 > snp_pos:
                    cur_pos += 1
                    if cur_pos >= len(l):
                        break
                    else:
                        snp_pos = int(l[cur_pos])
                if n_1 == snp_pos or n == snp_pos:
                    outFile.write(line)
                else:
                    pass
outFile.close()