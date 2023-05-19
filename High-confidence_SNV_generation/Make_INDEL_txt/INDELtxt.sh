#!/bin/bash
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

#Call Indels from the .maf files.
for f in *.maf
do
  python callINDELsFromMAF.py $f $f.indels
  echo "Finished processing file: $f"
done

#Concatenate by chromosome.
for i in {1..22}
do
  cat Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.*.indels > Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.ALL.maf.indels
  echo "Finished processing chromosome $i."
done

for i in {X,Y,MT}
do
  cat Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.*.indels > Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.ALL.maf.indels
  echo "Finished processing chromosome $i."
done

#Concatenate all the indels together
cat *ALL.maf.indels > all_indels_mixed.txt
echo "Merged all indels"

#Split the file for memory concerns
split -l 5000000 all_indels_mixed.txt indels
echo "Split Indels"

#Reformat for downstream filtering.  It is not really a bed format.  
for f in indelsa*
do
  python make_indel_bed.py $f
  echo "Finished making bedfile for: $f"
done
echo "Made bed files"

#Merge again
cat indelsa*chimp* > chimp_referenced_chp_hum_indels.bed
cat indelsa*human* > human_referenced_chp_hum_indels.bed
echo "Merged bed files"

#If commented out, need to remove _over from the +9 line.
sed 's/^/chr/' chimp_referenced_chp_hum_indels.bed > chimp_referenced_chp_hum_indels.lift
echo "Lifting Over"
./liftOver minMatch=0.9 chimp_referenced_chp_hum_indels.lift panTro4ToPanTro5.over.chain.gz chimp_referenced_chp_hum_indels_lifted.bed chimp_referenced_chp_hum_indels_lifted.err
sed 's/...//' chimp_referenced_chp_hum_indels_lifted.bed > chimp_referenced_chp_hum_indels_over.bed
echo "Done with liftover"

#Sort the files
sort human_referenced_chp_hum_indels.bed > all_indels_human_sorted.bed
sort chimp_referenced_chp_hum_indels_over.bed > all_indels_chimp_sorted.bed
echo "Sorted files"

#Filter out non-unique lines
uniq all_indels_human_sorted.bed > all_indels_human_sorted_uniq.bed
uniq all_indels_chimp_sorted.bed > all_indels_chimp_sorted_uniq.bed
echo "First pass filtering done"

#Switch the chromosome names
python switchChrom.py all_indels_chimp_sorted_uniq.bed all_indels_chimp_sorted_uniq_NC.bed chimpChroms.txt refseq
python switchChrom.py all_indels_human_sorted_uniq.bed all_indels_human_sorted_uniq_NC.bed humanChroms.txt refseq
echo "Chromosome names swapped"

#Remove duplicate indels
python rmdup_indels_hum.py
python rmdup_indels_nhp.py
echo "Second pass filtering done"

#Sorted again
sort all_indels_chimp_sorted_uniq_NC_rmdup.bed > all_indels_chimp_sorted_uniq_NC_rmdup_sorted_again.bed
sort all_indels_human_sorted_uniq_NC_rmdup.bed > all_indels_human_sorted_uniq_NC_rmdup_sorted_again.bed

#Reformat for WASP
python splitByChrom.py all_indels_human_sorted_uniq_rmdup_sorted_again.bed human indels huch
python splitByChrom.py all_indels_chimp_sorted_uniq_rmdup_sorted_again.bed chimp indels huch
echo "Reformatted for WASP, done"
