import sys
import os
import pandas as pd

f1 = sys.argv[1]
f2 = sys.argv[2]
f3 = sys.argv[3]
f4 = sys.argv[4]
r = sys.argv[5]
list_of_files = [f1, f2, f3, f4]
d_p = {}
out = open("All_Peaks_Merged_" + r + "_Fast.bed", 'w')
for f in list_of_files:
    o = open(f)
    d = {}
    for line in o:
        l = line.replace("\n", "").split("\t")
        if l[3] in d.keys():
            d[l[3]].append(l[12:])
            try:
                d_p[l[15].split("_")[1]].append(l[15])
            except:
                d_p[l[15].split("_")[1]] = [l[15]]
        else:
            d[l[3]] = [l[0:12], l[12:]]
            try:
                d_p[l[15].split("_")[1]].append(l[15])
            except:
                d_p[l[15].split("_")[1]] = [l[15]]
            try:
                d_p[l[3].split("_")[1]].append(l[3])
            except:
                d_p[l[3].split("_")[1]] = [l[3]]
    for key in d.keys():
        zipped = list(zip(*d[key]))
        new_peak = []
        for i in range(len(zipped)):
            #Add the chromsome, we know this is identical
            if i == 0:
                new_peak.append(zipped[i][0])
            #Add the minimum left value for the peak boundary
            elif i == 1:
                new_peak.append(min([float(x) for x in zipped[i]]))
            #Add the maximum right value for the peak boundary
            elif i == 2:
                new_peak.append(max([float(x) for x in zipped[i]]))
            #Lop off the annotations, we will add those back with the new peaklist 
            elif i == 3:
                cts = []
                jj = zipped[i][0].split("_")
                for j in zipped[i]:
                    cts.append(j.split("_")[1])
                new_peak_name = jj[0] + "_" + "_".join(jj[2:5]) + "_" + "_".join(cts)
                new_peak.append(new_peak_name)
            elif i in [4, 6, 7, 8, 9, 10]:
                new_peak.append(max([float(x) for x in zipped[i]]))
            elif i == 5:
                new_peak.append(".")
        out.write("\t".join([str(x) for x in new_peak])+"\n")
    o.close()

files = os.listdir()
for file in files:
    if "Peaklist_Final_Anno_" + r + ".bed" in file:
        v = pd.read_csv(file, sep = "\t", header = None)
        print(v)
        ll = d_p[file.split("_")[0]]
        v = v[~v[3].isin(ll)]
        for index, row in v.iterrows():
            l = list(row)
            lll = l[3].split("_")
            l[3] = lll[0] + "_" + "_".join(lll[2:5]) + "_" + lll[1]
            out.write("\t".join([str(x) for x in l]) + "\n")
out.close()

