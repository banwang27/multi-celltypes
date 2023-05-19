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

print("Merging overlapping peaks; be sure that the file is sorted in ascending order by chromosome and then the third column (indexed from 1)")
for index, row in z.iterrows():
    if len(prev):
        end = prev[2]
        beg = prev[1]
        if row[1] <= end and row[0] == prev[0]:
            prev[2] = row[2]
            ps.append(row[7])
            qs.append(row[8])
            ss.append(row[6])
            add_merge.append(row[3])
        else:
            prev[8] = max(qs)
            prev[7] = max(ps)
            prev[6] = max(ss)
            merged_peaks.append(add_merge)
            o.append(prev)
            prev = list(row)
            ss = [row[6]]
            ps = [row[7]]
            qs = [row[8]]
            add_merge = [row[3]]
    else:
        prev = list(row)
        ss = [row[6]]
        ps = [row[7]]
        qs = [row[8]]
        add_merge = [row[3]]
o.append(prev)
df = pd.DataFrame(o)
df.to_csv(a.replace(".narrowPeak", "_merged.narrowPeak").replace(".bed", "_merged.bed"), sep = "\t", index = False, header = False)
df2 = pd.DataFrame(merged_peaks)
df2.to_csv(a.replace(".narrowPeak", "_merged.csv").replace(".bed", "_merged.csv"), sep = "\t", index = False, header = False)
