import pandas as pd
import numpy as np
import sys

#closest is the file from bedtools closest, enh is the file of enhancer peaks, out_name is the prefix of the out file.
closest = sys.argv[1]
enh = sys.argv[2]
out_name = sys.argv[3]
vvv = pd.read_csv(closest, sep = "\t", header = None)

prev_enh = 0
d = {}
ind = 1
prev_gene = 0
for index, row in vvv.iterrows():
    #If it is the first time we are seeing this enhancer, we know that the gene is the closest and unique
    if row[3] != prev_enh:
        prev_enh = row[3]
        prev_gene = row[14]
        dist = abs((row[1] + row[2])/2 - (row[12] + row[13])/2)
        d[row[3]] = row[3] + "_" + row[14] +"_" + str(dist)
        ind = 1
    #Otherwise we need to be sure that the second closest gene is unique
    elif ind and row[13] != prev_gene:
        ind = 0
        dist = abs((row[1] + row[2])/2 - (row[12] + row[13])/2)
        d[row[3]] = d[row[3]] + "_" + row[14] +"_" + str(dist)
        
vvv = pd.read_csv(enh, sep = "\t", header = None)
out = []
for index, row in vvv.iterrows():
    l = list(row)
    l[3] = d[l[3]]
    out.append(l)
df = pd.DataFrame(out)
df.to_csv(out_name + "_Final_Enhancer.bed", sep = "\t", index = False, header = False)