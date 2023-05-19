#!/bin/bash
#SBATCH --time 48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

#From Encode pipeline, filter out based on replicates
#Need to do for humreffed and chpreffed
echo Filtering based on replicates
bedtools intersect -a ATAC7_CM_Hyb1_chpreffed.merged_peaks.narrowPeak -b ATAC7_CM_Hyb1_chpreffed.rmdup.remap.kept.merged.filt_peaks.narrowPeak -wo | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k3n -k2n | uniq | 
bedtools intersect -a stdin -b ATAC8_CM_Hyb2_chpreffed.rmdup.remap.kept.merged.filt_peaks.narrowPeak -wo | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k3n -k2n | uniq > CM_Pooled_Chpreffed.narrowPeak

bedtools intersect -a ATAC7_CM_Hyb1_humreffed.merged_peaks.narrowPeak -b ATAC7_CM_Hyb1_humreffed.rmdup.remap.kept.merged.filt_peaks.narrowPeak -wo | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k3n -k2n | uniq | 
bedtools intersect -a stdin -b ATAC8_CM_Hyb2_humreffed.rmdup.remap.kept.merged.filt_peaks.narrowPeak -wo | 
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k3n -k2n | uniq > CM_Pooled_Humreffed.narrowPeak

#Or if there are 4 replicates
#bedtools intersect -a stdin -b Rep2.narrowPeak -wo | 
#awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k3n -k2n | uniq |
#bedtools intersect -a stdin -b Rep3.narrowPeak -wo | 
#awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k3n -k2n | uniq |
#bedtools intersect -a stdin -b Rep4.narrowPeak -wo | 
#awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k3n -k2n | uniq > PooledInRep1AndRep2.narrowPeak

#Need to merge to eliminate redundant peaks
python merge.py CM_Pooled_Chpreffed.narrowPeak
python merge.py CM_Pooled_Humreffed.narrowPeak

#Need to do liftover.
./liftOver -bedPlus=3 -minMatch=0.1 CM_Pooled_Chpreffed_merged.narrowPeak panTro6ToHg38.over.chain.gz CM_Pooled_Chpreffed_Lifted.narrowPeak CM_Pooled_Chpreffed_Lifted.err
./liftOver -bedPlus=3 -minMatch=0.1 CM_Pooled_Humreffed_merged.narrowPeak hg38ToPanTro6.over.chain.gz CM_Pooled_Humreffed_Lifted.narrowPeak CM_Pooled_Humreffed_Lifted.err
./liftOver -bedPlus=3 -minMatch=0.1 CM_Pooled_Humreffed_Lifted.narrowPeak panTro6ToHg38.over.chain.gz CM_Pooled_Humreffed_Lifteded.narrowPeak CM_Pooled_Humreffed_Lifteded.err

echo Intersecting the lifted files
bedtools intersect -a CM_Pooled_Humreffed_Lifteded.narrowPeak -b CM_Pooled_Chpreffed_Lifted.narrowPeak -wo > CM_Peaklist_Intersect.narrowPeak

#Need to find intersecting peaks
python find_partners.py CM_Peaklist_Intersect.narrowPeak 0.25

#Filter out the blacklisted regions
#From Encode pipeline, removes blacklisted peaks
bedtools intersect -v -a CM_Peaklist_Intersect_filtered_union.narrowPeak -b hg38.blacklist.bed | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | sort -k1,1 -k3n -k2n > CM_Peaklist_Intersect_Union_Filt.narrowPeak
#Merge any overlapping peaks

python merge.py CM_Peaklist_Intersect_Union_Filt.narrowPeak

mv CM_Peaklist_Intersect_Union_Filt_merged.narrowPeak CM_Peaklist_Intersect_Union_Filt_Humreffed_merged.narrowPeak
#Lift back over to the chimpanzee reference
./liftOver -bedPlus=3 -minMatch=0.1 CM_Peaklist_Intersect_Union_Filt_Humreffed_merged.narrowPeak hg38ToPanTro6.over.chain.gz CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged.narrowPeak CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged.err

sort -k1,1 -k3n -k2n CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged.narrowPeak > CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted.narrowPeak

#Instead of merging need to merge and record any peaks that were merged.
#Output the list of merged peaks and filter those out of the human referenced and chimpanzee referenced files using the partners from find_partners
python merge.py CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted.narrowPeak

#Finally remove peaks that are not 1-1 orthologous between species
python filter_11.py CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted_merged.csv CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged.err CM_Peaklist_Intersect_Union_Filt_Humreffed_merged.narrowPeak CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted_merged.narrowPeak

#Now intersect with the list of promoters we have
bedtools intersect -a CM_Peaklist_Intersect_Union_Filt_Humreffed_merged_Final.narrowPeak -b Human_Promoters_Ortho_Sorted_hg38.bed -wa -wb > CM_Peaklist_Humreffed_Promoter_Intersect.bed
bedtools intersect -a CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted_merged_Final.narrowPeak -b Human_Promoters_Ortho_Sorted_PanTro6.bed -wa -wb > CM_Peaklist_Chpreffed_Promoter_Intersect.bed

#Next use a python script to do some annotation and split into promoters and putative enhancers
python anno_prom.py CM_Peaklist_Humreffed_Promoter_Intersect.bed CM_Peaklist_Intersect_Union_Filt_Humreffed_merged_Final.narrowPeak CM_Peaklist_Humreffed
python anno_prom.py CM_Peaklist_Chpreffed_Promoter_Intersect.bed CM_Peaklist_Intersect_Union_Filt_Chpreffed_merged_sorted_merged_Final.narrowPeak CM_Peaklist_Chpreffed

#Next use the sorted Orthologous TSS and bedtools closest to assign enhancers to genes (coupled to a python script)
#Can sort things appropriately with sort -k1,1 -k2,2n
sort -k1,1 -k2,2n CM_Peaklist_Humreffed_Temp_Enhancer.bed > CM_Peaklist_Humreffed_Temp_Enhancer_Sorted.bed
bedtools closest -k 10 -mdb all -a CM_Peaklist_Humreffed_Temp_Enhancer_Sorted.bed -b Human_Promoters_Ortho_Sorted_hg38.bed > CM_Peaklist_Temp_Enhancer_Closest_Humreffed.bed
sort -k1,1 -k2,2n CM_Peaklist_Chpreffed_Temp_Enhancer.bed > CM_Peaklist_Chpreffed_Temp_Enhancer_Sorted.bed
bedtools closest -k 10 -mdb all -a CM_Peaklist_Chpreffed_Temp_Enhancer_Sorted.bed -b Human_Promoters_Ortho_Sorted_PanTro6.bed > CM_Peaklist_Temp_Enhancer_Closest_Chpreffed.bed

#Another python script to annotate putative enhancers.
python anno_enh.py CM_Peaklist_Temp_Enhancer_Closest_Humreffed.bed CM_Peaklist_Humreffed_Temp_Enhancer_Sorted.bed CM_Peaklist_Humreffed
python anno_enh.py CM_Peaklist_Temp_Enhancer_Closest_Chpreffed.bed CM_Peaklist_Chpreffed_Temp_Enhancer_Sorted.bed CM_Peaklist_Chpreffed

#Recombine the two things
cat CM_Peaklist_Final_Enhancer_Humreffed.bed CM_Peaklist_Humreffed_Promoter_Final.bed > CM_Peaklist_Final_Anno_Humreffed.bed
cat CM_Peaklist_Chpreffed_Final_Enhancer.bed CM_Peaklist_Chpreffed_Final_Promoter.bed > CM_Peaklist_Final_Anno_Chpreffed.bed

#Convert to a gtf file
python np_to_gtf.py CM_Peaklist_Final_Anno_Humreffed.bed
python np_to_gtf.py CM_Peaklist_Final_Anno_Chpreffed.bed
