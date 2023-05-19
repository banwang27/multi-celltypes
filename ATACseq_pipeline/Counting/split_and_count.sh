#!/bin/bash
#SBATCH --time 48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

#Split the files, requires ASE_SNPs.FILTER.SPLIT_SPECIES_HUMAN.bed
for bamfile in *humreffed*.bam;
do
    echo splitting $bamfile
    samtools sort -n -o ${bamfile::-4}.sortn.bam $bamfile
    python splitSpeciesReads.py ASE_SNPs.FILTER.SPLIT_SPECIES_HUMAN.bed ${bamfile::-4}.sortn.bam ${bamfile::-4}.Hume.bam ${bamfile::-4}.Chpe.bam
    samtools sort -o ${bamfile::-4}.Hume.sort.bam ${bamfile::-4}.Hume.bam
    samtools index ${bamfile::-4}.Hume.sort.bam
    samtools sort -o ${bamfile::-4}.Chpe.sort.bam ${bamfile::-4}.Chpe.bam
    samtools index ${bamfile::-4}.Chpe.sort.bam
done

#Split the files, requires ASE_SNPs.FILTER.SPLIT_SPECIES_HUMAN.bed
for bamfile in *chpreffed*.bam;
do
    echo splitting $bamfile
    samtools sort -n -o ${bamfile::-4}.sortn.bam $bamfile
    python splitSpeciesReads.py ASE_SNPs.FILTER.SPLIT_SPECIES_CHIMP.bed ${bamfile::-4}.sortn.bam ${bamfile::-4}.Hume.bam ${bamfile::-4}.Chpe.bam
    samtools sort -o ${bamfile::-4}.Hume.sort.bam ${bamfile::-4}.Hume.bam
    samtools index ${bamfile::-4}.Hume.sort.bam
    samtools sort -o ${bamfile::-4}.Chpe.sort.bam ${bamfile::-4}.Chpe.bam
    samtools index ${bamfile::-4}.Chpe.sort.bam
done

#Counts reads in peaks with HTSeq
python -m HTSeq.scripts.count -t peak -i peak_ID -s no -m union -r pos ATAC7_CM_Hyb1_humreffed.rmdup.remap.kept.merged.filt.Chpe.sort.bam CM_Peaklist_Final_Anno_Humreffed.gtf > ATAC7_CM_Hyb1_Chimp_humreffed.txt
python -m HTSeq.scripts.count -t peak -i peak_ID -s no -m union -r pos ATAC7_CM_Hyb1_humreffed.rmdup.remap.kept.merged.filt.Hume.sort.bam CM_Peaklist_Final_Anno_Humreffed.gtf > ATAC7_CM_Hyb1_Human_humreffed.txt
python -m HTSeq.scripts.count -t peak -i peak_ID -s no -m union -r pos ATAC8_CM_Hyb2_humreffed.rmdup.remap.kept.merged.filt.Chpe.sort.bam CM_Peaklist_Final_Anno_Humreffed.gtf > ATAC8_CM_Hyb2_Chimp_humreffed.txt
python -m HTSeq.scripts.count -t peak -i peak_ID -s no -m union -r pos ATAC8_CM_Hyb2_humreffed.rmdup.remap.kept.merged.filt.Hume.sort.bam CM_Peaklist_Final_Anno_Humreffed.gtf > ATAC8_CM_Hyb2_Human_humreffed.txt

#Might have to delete non traditional chromosomes (e.g. chrM, chrUn) although they should be filtered out earlier.
python -m HTSeq.scripts.count -t peak -i peak_ID -s no -m union -r pos ATAC7_CM_Hyb1_chpreffed.rmdup.remap.kept.merged.filt.Chpe.sort.bam CM_Peaklist_Final_Anno_Chpreffed.gtf > ATAC7_CM_Hyb1_Chimp_chpreffed.txt
python -m HTSeq.scripts.count -t peak -i peak_ID -s no -m union -r pos ATAC7_CM_Hyb1_chpreffed.rmdup.remap.kept.merged.filt.Hume.sort.bam CM_Peaklist_Final_Anno_Chpreffed.gtf > ATAC7_CM_Hyb1_Human_chpreffed.txt
python -m HTSeq.scripts.count -t peak -i peak_ID -s no -m union -r pos ATAC8_CM_Hyb2_chpreffed.rmdup.remap.kept.merged.filt.Chpe.sort.bam CM_Peaklist_Final_Anno_Chpreffed.gtf > ATAC8_CM_Hyb2_Chimp_chpreffed.txt
python -m HTSeq.scripts.count -t peak -i peak_ID -s no -m union -r pos ATAC8_CM_Hyb2_chpreffed.rmdup.remap.kept.merged.filt.Hume.sort.bam CM_Peaklist_Final_Anno_Chpreffed.gtf > ATAC8_CM_Hyb2_Human_chpreffed.txt
