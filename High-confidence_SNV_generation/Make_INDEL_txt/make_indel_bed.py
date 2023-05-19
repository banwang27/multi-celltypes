import os
import pandas as pd

def complement(bp):
    if bp == "A":
        return("T")
    elif bp == "T":
        return("A")
    elif bp == "C":
        return("G")
    elif bp == "-":
        return "-"
    else:
        return("C")


out_chp = []
out_hum = []
csv = pd.read_csv("all_indels_mixed.txt", sep = "\t")
for index, row in csv.iterrows():
    new_row_hum = [row[0], row[1]-1, row[1], row[3]+"|"+row[4]]
    if row[2] == row[5]:
    	new_row_chp = [row[7], row[6]-1, row[6], row[3]+"|"+row[4]]
    else:
    	new_row_chp = [row[7], row[6]-1, row[6], complement(row[3])+"|"+complement(row[4])]
    out_hum.append(new_row_hum)
    out_chp.append(new_row_chp)
df_hum = pd.DataFrame(out_hum)
df_hum.to_csv("human_referenced_chp_hum_indels.bed", sep = "\t", header = False, index = False)
df_chp = pd.DataFrame(out_chp)
df_chp.to_csv("chimp_referenced_chp_hum_indels.bed", sep = "\t", header = False, index = False)