import pandas as pd
import numpy as np
import sys

#prom_inter is the file from bedtools intersect, peaklist is the list of peaks, and out is the prefix of the output file
prom_inter = sys.argv[1]
peaklist = sys.argv[2]
out = sys.argv[3]
v = pd.read_csv(peaklist, sep = "\t", header = None)
v2 = pd.read_csv(prom_inter, sep = "\t", header = None)

#Preprocess to remove duplicate promoters
d = {}
keep = []
proms = []
for index, r in v2.iterrows():
    prom_id = r[11] + str(r[12]) + str(r[13]) + str(r[14])
    proms.append(r[3])
    if prom_id in d.keys():
        d[prom_id].append(list(r))
    else:
        d[prom_id] = [list(r)]

#If the combined peak boundaries are wider than the promoter, we keep the wider interval
#Basically sticks with the whole union strategy
for key in d.keys():
    g = d[key]
    if len(d[key]) == 1:
        i = g[0]
        keep.append(i[0:12] + [min([i[1], i[12]]), max([i[2], i[13]])] + [i[14]])
    else:
        ss = []
        ps = []
        qs = []
        left = []
        right = []
        for i in g:
            ss.append(i[6])
            ps.append(i[7])
            qs.append(i[8])
            left = left + [i[1], i[12]]
            right = right + [i[2], i[13]]
        keep.append(i[0:6] + [max(ss), max(ps), max(qs)] + i[9:12] + [min(left), max(right)] + [i[14]])
v2 = pd.DataFrame(keep)

v2 = v2.set_index(3)
v = v.set_index(3)

#Need to do something if the promoters overlap
#So we just assign the peak to both promoters
d = {}
keep2 = []
for index, row in v2.iterrows():
    if index not in d.keys():
        d[index] = [list(row)]
    else:
        d[index].append(list(row))

for key in d.keys():
    if len(d[key]) < 2:
        keep2.append(d[key][0][0:3] + [key] + d[key][0][3:])
    else:
        genes = []
        lefts = []
        rights = []
        for i in d[key]:
            gene = i[13]
            if gene not in genes:
                genes.append(gene)
                lefts.append(i[11])
                rights.append(i[12])
        new_add = i[0:3] + [key] + i[3:11] + [min(lefts), max(rights)] + ["_".join(genes)]
        keep2.append(new_add)
v2 = pd.DataFrame(keep2)
v2 = v2.set_index(3)

#Go through and replace peaks that overlap promoters with the unionized promoter to maximize statistical power
out_prom = []
out_enh = []
for index, row in v.iterrows():
    if index in list(v2.index):
        r = v2.loc[index]
        prom_id = r[11] + str(r[12]) + str(r[13]) + str(r[14])
        out_prom.append(list([r[11], r[12], r[13]] + [index + "_promoter_" + r[14]] + list(row)[3:5] + list(r[5:8]) + list(row)[8:]))
    else:
        if index not in proms:
            out_enh.append(list(row)[0:3] + [index + "_enhancer"] + list(row)[3:])
            
df = pd.DataFrame(out_prom)
df2 = pd.DataFrame(out_enh)
df2.to_csv(out + "_Temp_Enhancer.bed", sep = "\t", header = False, index = False)
df.to_csv(out + "_Final_Promoter.bed", sep = "\t", header = False, index = False)

