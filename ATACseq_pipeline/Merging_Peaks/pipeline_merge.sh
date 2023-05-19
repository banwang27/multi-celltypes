#!/bin/bash
#SBATCH --time 48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

#First, we iteratively intersect all files removing the -a file as we go.
#This will introduce redundancy
bedtools intersect -a CM_Peaklist_Final_Anno_Humreffed.bed -b HP_Peaklist_Final_Anno_Humreffed.bed MN_Peaklist_Final_Anno_Humreffed.bed PP_Peaklist_Final_Anno_Humreffed.bed SKM_Peaklist_Final_Anno_Humreffed.bed -wo > All_Intersect_Peaklist.bed
bedtools intersect -a HP_Peaklist_Final_Anno_Humreffed.bed -b MN_Peaklist_Final_Anno_Humreffed.bed PP_Peaklist_Final_Anno_Humreffed.bed SKM_Peaklist_Final_Anno_Humreffed.bed -wo > All_Intersect_Peaklist_2.bed
bedtools intersect -a MN_Peaklist_Final_Anno_Humreffed.bed -b PP_Peaklist_Final_Anno_Humreffed.bed SKM_Peaklist_Final_Anno_Humreffed.bed -wo > All_Intersect_Peaklist_3.bed
bedtools intersect -a PP_Peaklist_Final_Anno_Humreffed.bed -b SKM_Peaklist_Final_Anno_Humreffed.bed -wo > All_Intersect_Peaklist_4.bed

#This script adds all the peaks together and adds back any peaks that lacked a partner in any of the above files
#You may need to add a column in front of the second columns of chromsome names to prevent errors in column acccess
python deal_fast.py All_Intersect_Peaklist.bed All_Intersect_Peaklist_2.bed All_Intersect_Peaklist_3.bed All_Intersect_Peaklist_4.bed Humreffed

#Repeat for chpreffed
bedtools intersect -a CM_Peaklist_Final_Anno_Chpreffed.bed -b HP_Peaklist_Final_Anno_Chpreffed.bed MN_Peaklist_Final_Anno_Chpreffed.bed PP_Peaklist_Final_Anno_Chpreffed.bed SKM_Peaklist_Final_Anno_Chpreffed.bed -wo > All_Intersect_Peaklist_Chpreffed.bed
bedtools intersect -a HP_Peaklist_Final_Anno_Chpreffed.bed -b MN_Peaklist_Final_Anno_Chpreffed.bed PP_Peaklist_Final_Anno_Chpreffed.bed SKM_Peaklist_Final_Anno_Chpreffed.bed -wo > All_Intersect_Peaklist_Chpreffed_2.bed
bedtools intersect -a MN_Peaklist_Final_Anno_Chpreffed.bed -b PP_Peaklist_Final_Anno_Chpreffed.bed SKM_Peaklist_Final_Anno_Chpreffed.bed -wo > All_Intersect_Peaklist_Chpreffed_3.bed
bedtools intersect -a PP_Peaklist_Final_Anno_Chpreffed.bed -b SKM_Peaklist_Final_Anno_Chpreffed.bed -wo > All_Intersect_Peaklist_Chpreffed_4.bed

python deal_fast.py All_Intersect_Peaklist_Chpreffed.bed All_Intersect_Peaklist_Chpreffed_2.bed All_Intersect_Peaklist_Chpreffed_3.bed All_Intersect_Peaklist_Chpreffed_4.bed Chpreffed

#Sort
sort -k1,1 -k3n -k2n All_Peaks_Merged_Humreffed_Fast.bed > All_Peaks_Merged_Humreffed_Sorted.narrowPeak
sort -k1,1 -k3n -k2n All_Peaks_Merged_Chpreffed_Fast.bed > All_Peaks_Merged_Chpreffed_Sorted.narrowPeak

#Need to merge to eliminate redundant peaks
#This is a special merge operation that instead seeks to include all the cell types with an intersecting peak
#It also adds a unique tag to the peak to prevent peak name intersection
python merge_all.py All_Peaks_Merged_Humreffed_Sorted.narrowPeak
python merge_all.py All_Peaks_Merged_Chpreffed_Sorted.narrowPeak

#Need to do reciprocal liftover to eliminate peaks without clear orthologs.
./liftOver -bedPlus=3 -minMatch=0.1 All_Peaks_Merged_Chpreffed_Sorted_merged.narrowPeak panTro6ToHg38.over.chain.gz All_Peaks_Chpreffed_Lifted.narrowPeak All_Peaks_Chpreffed_Lifted.err
./liftOver -bedPlus=3 -minMatch=0.1 All_Peaks_Merged_Humreffed_Sorted_merged.narrowPeak hg38ToPanTro6.over.chain.gz All_Peaks_Humreffed_Lifted.narrowPeak All_Peaks_Humreffed_Lifted.err
./liftOver -bedPlus=3 -minMatch=0.1 All_Peaks_Humreffed_Lifted.narrowPeak panTro6ToHg38.over.chain.gz All_Peaks_Humreffed_Lifteded.narrowPeak All_Peaks_Humreffed_Lifteded.err

#echo Intersecting the lifted files
bedtools intersect -a All_Peaks_Humreffed_Lifteded.narrowPeak -b All_Peaks_Chpreffed_Lifted.narrowPeak -wo > All_Peaks_Peaklist_Intersect.narrowPeak

#Need to find intersecting peaks
python find_partners_all.py All_Peaks_Peaklist_Intersect.narrowPeak 0.25

#Filter out the blacklisted regions
#From Encode pipeline, removes blacklisted peaks
#https://www.nature.com/articles/s41598-019-45839-z for blacklist details
bedtools intersect -v -a All_Peaks_Peaklist_Intersect_filtered_union.narrowPeak -b hg38.blacklist.bed | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | sort -k1,1 -k3n -k2n > All_Peaks_Peaklist_Intersect_Union_Filt.narrowPeak
#Merge any overlapping peaks
sort -k1,1 -k3n -k2n All_Peaks_Peaklist_Intersect_Union_Filt.narrowPeak > All_Peaks_Peaklist_Intersect_Union_Filt_sorted.narrowPeak
python merge_all2.py All_Peaks_Peaklist_Intersect_Union_Filt_sorted.narrowPeak

mv All_Peaks_Peaklist_Intersect_Union_Filt_merged.narrowPeak All_Peaks_Peaklist_Intersect_Union_Filt_Humreffed_merged.narrowPeak
#Lift back over to the chimpanzee reference
./liftOver -bedPlus=3 -minMatch=0.1 All_Peaks_Peaklist_Intersect_Union_Filt_Humreffed_sorted_merged.narrowPeak hg38ToPanTro6.over.chain.gz All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged.narrowPeak All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged.err

sort -k1,1 -k3n -k2n All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged.narrowPeak > All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted.narrowPeak
#Instead of merging need to merge and record any peaks that were merged.
#Output the list of merged peaks and filter those out of the human referenced and chimpanzee referenced files using the partners from find_partners
python merge_all2.py All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted.narrowPeak

#Finally remove peaks that are not 1-1 orthologous between species after liftover
#This can occur due to liftover leading to two peaks being merged, deleted, etc.
python filter_11.py All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted_merged.csv All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged.err All_Peaks_Peaklist_Intersect_Union_Filt_Humreffed_sorted_merged.narrowPeak All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted_merged.narrowPeak

#Now intersect with the list of promoters we have
bedtools intersect -a All_Peaks_Peaklist_Intersect_Union_Filt_Humreffed_sorted_merged_Final.narrowPeak -b Human_Promoters_Ortho_Sorted_hg38.bed -wa -wb > All_Peaks_Peaklist_Humreffed_Promoter_Intersect.bed
bedtools intersect -a All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted_merged_Final.narrowPeak -b Human_Promoters_Ortho_Sorted_PanTro6.bed -wa -wb > All_Peaks_Peaklist_Chpreffed_Promoter_Intersect.bed

#Next use a python script to do some annotation and split into promoters and putative enhancers
python anno_prom.py All_Peaks_Peaklist_Humreffed_Promoter_Intersect.bed All_Peaks_Peaklist_Intersect_Union_Filt_Humreffed_sorted_merged_Final.narrowPeak All_Peaks_Peaklist_Humreffed
python anno_prom.py All_Peaks_Peaklist_Chpreffed_Promoter_Intersect.bed All_Peaks_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted_merged_Final.narrowPeak All_Peaks_Peaklist_Chpreffed

#Next use the sorted Orthologous TSS and bedtools closest to assign enhancers to genes (coupled to a python script)
#Can sort things appropriately with sort -k1,1 -k2,2n
sort -k1,1 -k2,2n All_Peaks_Peaklist_Humreffed_Temp_Enhancer.bed > All_Peaks_Peaklist_Humreffed_Temp_Enhancer_Sorted.bed
bedtools closest -k 10 -mdb all -a All_Peaks_Peaklist_Humreffed_Temp_Enhancer_Sorted.bed -b Human_Promoters_Ortho_Sorted_hg38.bed > All_Peaks_Peaklist_Temp_Enhancer_Closest_Humreffed.bed
sort -k1,1 -k2,2n All_Peaks_Peaklist_Chpreffed_Temp_Enhancer.bed > All_Peaks_Peaklist_Chpreffed_Temp_Enhancer_Sorted.bed
bedtools closest -k 10 -mdb all -a All_Peaks_Peaklist_Chpreffed_Temp_Enhancer_Sorted.bed -b Human_Promoters_Ortho_Sorted_PanTro6.bed > All_Peaks_Peaklist_Temp_Enhancer_Closest_Chpreffed.bed

#Another python script to annotate again.
python anno_enh.py All_Peaks_Peaklist_Temp_Enhancer_Closest_Humreffed.bed All_Peaks_Peaklist_Humreffed_Temp_Enhancer_Sorted.bed All_Peaks_Peaklist_Humreffed
python anno_enh.py All_Peaks_Peaklist_Temp_Enhancer_Closest_Chpreffed.bed All_Peaks_Peaklist_Chpreffed_Temp_Enhancer_Sorted.bed All_Peaks_Peaklist_Chpreffed

#Recombine the two things
cat All_Peaks_Peaklist_Humreffed_Final_Enhancer.bed All_Peaks_Peaklist_Humreffed_Final_Promoter.bed > All_Peaks_Peaklist_Final_Anno_Humreffed.bed
cat All_Peaks_Peaklist_Chpreffed_Final_Enhancer.bed All_Peaks_Peaklist_Chpreffed_Final_Promoter.bed > All_Peaks_Peaklist_Final_Anno_Chpreffed.bed

#Convert to a gtf file
python np_to_gtf.py All_Peaks_Peaklist_Final_Anno_Humreffed.bed
python np_to_gtf.py All_Peaks_Peaklist_Final_Anno_Chpreffed.bed
