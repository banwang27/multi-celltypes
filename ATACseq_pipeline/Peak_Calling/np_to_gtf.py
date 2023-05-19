import sys
import pandas as pd
import numpy as np

#convert a narrow peak format file to a gtf file
f = sys.argv[1]
file = open(f)
out = open(f.replace(".bed", ".gtf").replace(".narrowPeak", ".gtf"), 'w')
for line in file:
    l = line.split("\t")
    to_write = [l[0], "macs2", "peak", l[1], l[2], l[4], l[5], "."]
    ninth_field = "peak_ID " + '\"' + l[3] + '\"' + ";" + "p-value " + '\"' + l[7] + '\"' + ";" + "q-value " + '\"' + l[8] + '\"' + ";\n"
    to_write = to_write + [ninth_field]
    out.write("\t".join(to_write))
file.close()
out.close()