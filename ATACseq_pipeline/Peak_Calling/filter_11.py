import sys
import pandas as pd
import numpy as np

#Merge is the csv file from the final merge step, err is the error file from liftover, humreffed and chpreffed are the peaks from the respective sides
merge = sys.argv[1]
err = sys.argv[2]
humreffed = sys.argv[3]
chpreffed = sys.argv[4]

#Remove things that became intersecting after final liftover
bad = []
for index, row in pd.read_csv(merge, sep = "\t", header = None).iterrows():
    l = list(row)
    if len(l) > 2:
        for i in l:
            bad.append(i)
#Remove things that failed the final liftover
f = open(err)
for line in f:
    if "#" not in line:
        l = line.split("\t")[3]
        bad.append(l)
f.close()

#Actually remove those things from humreffed
out = open(humreffed.replace(".narrowPeak", "_Final.narrowPeak"), 'w')
f = open(humreffed)
for line in f:
    l = line.split("\t")
    if l[3] not in bad:
        out.write(line)
out.close()
f.close()

#Actually remove those things from chpreffed
out = open(chpreffed.replace(".narrowPeak", "_Final.narrowPeak").replace(".bed", "_Final.bed"), 'w')
f = open(chpreffed)
for line in f:
    l = line.split("\t")
    if l[3] not in bad:
        out.write(line)
out.close()
f.close()
