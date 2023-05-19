#!/bin/bash
#SBATCH --time 48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G


for file in *chpreffed*.bam
do
    samtools view -b $file chr1 chr2A chr2B chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${file::-4}.filt.bam
done

for file in *humreffed*.bam
do
    samtools view -b $file chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${file::-4}.filt.bam
done

#Merge, can add more arguments if > 2 replicates
samtools merge ATAC7_CM_Hyb1_chpreffed.merged.bam ATAC7_CM_Hyb1_chpreffed.rmdup.remap.kept.merged.filt.filt.bam ATAC8_CM_Hyb2_chpreffed.rmdup.remap.kept.merged.filt.filt.bam
samtools merge ATAC7_CM_Hyb1_humreffed.merged.bam ATAC7_CM_Hyb1_humreffed.rmdup.remap.kept.merged.filt.filt.bam ATAC8_CM_Hyb2_humreffed.rmdup.remap.kept.merged.filt.filt.bam

#We convert to a bed file as this is recommended by various people.
#Similar arguments to encode are used.
for file in *.bam
do
	bedtools bamtobed -i $file > ${file::-4}.bed
	macs2 callpeak -t ${file::-4}.bed -f BED -p 0.01 --nomodel --shift 75 --extsize 150 -B --SPMR --keep-dup all --call-summits -n ${file::-4}
done
