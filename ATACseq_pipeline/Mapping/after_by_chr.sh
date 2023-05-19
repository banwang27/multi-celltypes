#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

samtools merge rmdup.remap.kept.bam rmdup.chr1.kept.sort.bam rmdup.chr2A.kept.sort.bam rmdup.chr2B.kept.sort.bam rmdup.chr3.kept.sort.bam rmdup.chr4.kept.sort.bam rmdup.chr5.kept.sort.bam rmdup.chr6.kept.sort.bam rmdup.chr7.kept.sort.bam rmdup.chr8.kept.sort.bam rmdup.chr9.kept.sort.bam rmdup.chr10.kept.sort.bam rmdup.chr11.kept.sort.bam rmdup.chr12.kept.sort.bam rmdup.chr13.kept.sort.bam rmdup.chr14.kept.sort.bam rmdup.chr15.kept.sort.bam rmdup.chr16.kept.sort.bam rmdup.chr17.kept.sort.bam rmdup.chr18.kept.sort.bam rmdup.chr19.kept.sort.bam rmdup.chr20.kept.sort.bam rmdup.chr21.kept.sort.bam rmdup.chr22.kept.sort.bam rmdup.chrX.kept.sort.bam rmdup.chrY.kept.sort.bam
samtools sort rmdup.remap.kept.bam -o rmdup.remap.kept.sort.bam
samtools sort rmdup.keep.bam -o rmdup.keep.sort.bam
samtools merge rmdup.remap.kept.merged.bam rmdup.remap.kept.sort.bam rmdup.keep.sort.bam
samtools view -b -q 10 rmdup.remap.kept.merged.bam > rmdup.remap.kept.merged.filt.bam