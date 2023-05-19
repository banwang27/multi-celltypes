#!/bin/bash
#SBATCH --time 48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G

#ml load system rclone
#rclone copy fraser:backup/astarr/Myriad_Hybrid_ATAC_Final_Bams/ATAC16_HP_Hyb1_humreffed.rmdup.remap.kept.merged.filt.bam ./

#for file in *chpreffed*.bam
#do
    #samtools index $file
#    samtools view -b $file chr1 chr2A chr2B chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${file::-4}.filt.bam
#done

#for file in *humreffed*.bam
#do
    #samtools index $file
#    samtools view -b $file chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${file::-4}.filt.bam
#done

java -jar picard.jar SortSam I=ATAC11_SKM_Hyb1_chpreffed.rmdup.remap.kept.merged.filt.filt.bam O=ATAC11_SKM_Hyb1_chpreffed.rmdup.remap.kept.merged.filt.filt.sort.bam SORT_ORDER=queryname
java -jar picard.jar SortSam I=ATAC11_SKM_Hyb1_humreffed.rmdup.remap.kept.merged.filt.filt.bam O=ATAC11_SKM_Hyb1_humreffed.rmdup.remap.kept.merged.filt.filt.sort.bam SORT_ORDER=queryname

java -jar picard.jar SplitSamByNumberOfReads I=ATAC11_SKM_Hyb1_chpreffed.rmdup.remap.kept.merged.filt.filt.sort.bam OUTPUT=/scratch/users/astarr97/ATAC_HP SPLIT_TO_N_FILES=2 OUT_PREFIX=ATAC_SKM_SPLIT_CHPR
java -jar picard.jar SplitSamByNumberOfReads I=ATAC11_SKM_Hyb1_humreffed.rmdup.remap.kept.merged.filt.filt.sort.bam OUTPUT=/scratch/users/astarr97/ATAC_HP SPLIT_TO_N_FILES=2 OUT_PREFIX=ATAC_SKM_SPLIT_HUMR
