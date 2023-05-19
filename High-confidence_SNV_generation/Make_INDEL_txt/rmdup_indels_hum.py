import pandas as pd
csv = pd.read_csv("all_indels_human_sorted_uniq_NC.bed", sep = "\t", header = None)
output = []
l = []
for index, row in csv.iterrows():
    if not l:
        l.append(list(row))
    else:
        if l[0][1] == row[1]:
            l.append(list(row))
        else:
            i = len(l) - 1
            output.append(l[i])
            l = [list(row)]
df = pd.DataFrame(output)
df.to_csv("all_indels_human_sorted_uniq_NC_rmdup.bed", index = False, sep = "\t", header = False)