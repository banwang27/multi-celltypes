import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pandas as pd
import sys
import numpy as np

#Retains the lowest q-value associated with a peak (max log value).
a = sys.argv[1]
prev = []
qs = []
o = []
z = pd.read_csv(a, sep = "\t", header=None)
#z = z.sort_values([0,1])
merged_peaks = []
suffs = []
ss = []
ps = []
print("Merging overlapping peaks; be sure that the file is sorted in ascending order by chromosome and then the third column (indexed from 1)")
c = 0
for index, row in z.iterrows():
    if len(prev):
        c += 1
        end = prev[2]
        beg = prev[1]
        if row[1] <= end and row[0] == prev[0]:
            prev[2] = row[2]
            ps.append(row[7])
            qs.append(row[8])
            ss.append(row[6])
            suffs.append(row[3])
            add_merge.append(row[3])
        else:
            prev[8] = max(qs)
            prev[7] = max(ps)
            prev[6] = max(ss)
            tot = []
            for i in suffs:
                tot = tot + i.split("_")[4:]
            tot = list(set(tot))
            tot.sort()
            prev[3] = "peak" + "_" + "_".join(tot) + "_" + str(c)
            merged_peaks.append(add_merge)
            o.append(prev)
            prev = list(row)
            ss = [row[6]]
            ps = [row[7]]
            qs = [row[8]]
            suffs = [prev[3]]
            add_merge = [row[3]]
    else:
        prev = list(row)
        ss = [row[6]]
        ps = [row[7]]
        qs = [row[8]]
        add_merge = [row[3]]
        suffs = [prev[3]]

prev[8] = max(qs)
prev[7] = max(ps)
prev[6] = max(ss)
tot = []
for i in suffs:
    tot = tot + i.split("_")[4:]
tot = list(set(tot))
tot.sort()
prev[3] = "peak" + "_" + "_".join(tot) + "_" + str(c)
merged_peaks.append(add_merge)
o.append(prev)
prev = list(row)
ss = [row[6]]
ps = [row[7]]
qs = [row[8]]
suffs = [prev[3]]
add_merge = [row[3]]
        

df = pd.DataFrame(o)
df[1] = df[1].astype(int)
df[2] = df[2].astype(int)
df.to_csv(a.replace(".narrowPeak", "_merged.narrowPeak"), sep = "\t", index = False, header = False)
df2 = pd.DataFrame(merged_peaks)
df2.to_csv(a.replace(".narrowPeak", "_merged.csv"), sep = "\t", index = False, header = False)
