import os
import sys

directory = "."
species = sys.argv[1]
min_reads = sys.argv[2]
# All files in directory
allFiles = os.listdir(directory)

# Only the files we care about
files = []

for i in allFiles:
    if i.startswith("SNPCounts") and species in i:
        files.append(i)
snp_dict = {}
counter = 1
total = str(len(files))
d = {}
chrs = ["chr2A", "chr2B", "chrX", "chrY"]
for i in range(1, 23):
        chrs.append("chr" + str(i))
chrs = ["chr1"]
for i in files:
    with open(i, 'r') as inFile:
        ind = 0
        for line in inFile:
            if ind:
                chrom, pos, ref, alt, non = line.strip().split('\t')
                if chrom not in d.keys():
                    d[chrom] = set()
                ref, alt, non = [int(j) for j in [ref,alt,non]]
                if ref + alt >= int(min_reads):
                    d[chrom].add(int(pos))
            else:
                ind += 1
for chrom in d.keys():
    if chrom in chrs:
        l = list(d[chrom])
        l.sort()
        outFile = open(chrom + "_" + species + "MIN_" + min_reads + "_READS" + "_" + "KEEP.txt", "w+")
        for j in l:
            outFile.write(str(j))
            outFile.write("\n")
        outFile.close()