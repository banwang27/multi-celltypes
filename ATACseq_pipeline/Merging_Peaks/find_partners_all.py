import os
import pandas as pd
import sys
import numpy as np

#Read in the intersected narrowPeak file derived from liftover
#cutoff is the argument that determines the degree of overlap required between two peaks for them to be considered partners.
y = sys.argv[1]
cutoff = float(sys.argv[2])
m = pd.read_csv(sys.argv[1], header = None, sep = "\t")
print(m.head())

#Read in as dictionary as there can be multiple peaks that overlap
d1 = {}
for index, row in m.iterrows():
    k = str(row[1]) + "-" + str(row[2])
    if k in d1.keys():
        d1[k].append(list(row))
    else:
        d1[k] = [list(row)]
        
def union(p1, p2):
    return [min(p1[0], p2[0]), max(p1[1], p2[1])]

#We keep the human reference values for the peak signal, p, q value as well as the cell types in which the peak was called.
out = []
for key in d1.keys():
    cur = d1[key]
    #If the length is 1 and there is sufficient overlap, then we take the union of the intervals and keep the peak.
    if len(cur) == 1:
        z = cur[0]
        if z[20]/(z[2] - z[1]) > cutoff or z[20]/(z[12] - z[11]) > cutoff:
            u = union([z[1], z[2]], [z[11], z[12]])
            #print([z[1], z[2]], [z[11], z[12]])
            #print(u)
            z[12] = u[1]
            z[11] = u[0]
            #Because humreffed goes first!
            z[13] = z[3]
            out.append(z)
    else:
        temp = []
        for z in cur:
            if z[20]/(z[2] - z[1]) > cutoff or z[20]/(z[12] - z[11]) > cutoff:
                u = union([z[1], z[2]], [z[11], z[12]])
                #print([z[1], z[2]], [z[11], z[12]])
                #print(u)
                
                z[12] = u[1]
                z[11] = u[0]
                #Because humreffed goes first!
                z[16] = z[6]
                z[17] = z[7]
                z[18] = z[8]
                z[13] = z[3]
                temp.append(z)
        if temp:
            if len(temp) == 1:
                out.append(temp[0])
            else:
                mins = []
                maxs = []
                maxq = []
                maxp = []
                maxsig = []
                tot = []
                for j in temp:
                    mins.append(j[11])
                    maxs.append(j[12])
                    maxq.append(j[8])
                    maxp.append(j[7])
                    maxsig.append(j[6])
                    tot = tot + j[3].split("_")[1:len(j[3].split("_"))-1]
                j = temp[0]
                j[11] = min(mins) 
                j[12] = max(maxs)
                j[18] = max(maxq)
                j[17] = max(maxp)
                j[16] = max(maxsig)
                tot = list(set(tot))
                tot.sort()
                j[13] = j[3].split("_")[0] + "_".join(tot) + str(j[len(j[3].split("_"))-1])
                out.append(j)
out_real = []
partners = []
for o in out:
    out_real.append(list(o[10:]))
    partners.append([o[3], o[13]])
df = pd.DataFrame(out_real)
df.to_csv(y.replace(".narrowPeak", "_filtered_union.narrowPeak").replace(".bed", "_filtered_union.bed"), sep = "\t", index = False, header = False)
df = pd.DataFrame(partners)
df.to_csv(y.replace(".narrowPeak", "_partnered_union.narrowPeak").replace(".bed", "_partnered_union.bed"), sep = "\t", index = False, header = False)     