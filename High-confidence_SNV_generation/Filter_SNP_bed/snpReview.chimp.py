#!/usr/bin/env python

import os

directory = "."

# All files in directory
allFiles = os.listdir(directory)

# Only the files we care about
files = []

for i in allFiles:
    if i.startswith("snpCounts") and "chimp" in i:
        files.append(i)
print(files)
snp_dict = {}
counter = 1
total = str(len(files))

for i in files:
    with open(i, 'r') as inFile:
        inFile.readline()
        print("Reading file: " + i)
        for line in inFile:
            chrom, pos, ref, alt, non = line.strip().split('\t')

            ref, alt, non = [int(j) for j in [ref,alt,non]]

            if (ref + alt + non) == 0:
                continue  

            if "H2" in i:
                h_correct = ref
                h_incorrect = alt
                h_none = non
                c_correct = 0
                c_incorrect = 0
                c_none = 0

            elif "C3" in i:
                c_correct = alt
                c_incorrect = ref
                c_none = non
                h_correct = 0
                h_incorrect = 0
                h_none = 0

            else:
                print("Error in filename: " + i)
                break

            if chrom in snp_dict:
                if pos in snp_dict[chrom]:
                    snp_dict[chrom][pos] = [sum(x) for x in zip(snp_dict[chrom][pos], [h_correct, h_incorrect, h_none, c_correct, c_incorrect, c_none])]
                else:
                    snp_dict[chrom][pos] = [h_correct, h_incorrect, h_none, c_correct, c_incorrect, c_none]

            else:
                snp_dict[chrom] = {}
                snp_dict[chrom][pos] = [h_correct, h_incorrect, h_none, c_correct, c_incorrect, c_none]

    print("Finished file " + str(counter) + " of " + total)
    counter += 1

out = open("PanTro5_SNPreview.ALL.SPLIT_SPECIES.CNCC.txt", 'w')

print("Writing output...")
for i in sorted(snp_dict.keys()):
	for j in sorted(snp_dict[i].keys()):
		vals = [str(i) for i in snp_dict[i][j]]
		outlist = [i,j] + vals
		out.write('\t'.join(outlist) + '\n')


out.close()
print("Finished")
