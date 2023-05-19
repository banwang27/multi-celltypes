#!/bin/bash
#SBATCH --time 48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G


#We convert to a bed file as this is recommended by various people.
#Similar arguments to encode are used.
for file in *filt.filt.bam
do
	bedtools bamtobed -i $file > ${file::-4}.bed
	macs2 callpeak -t ${file::-4}.bed -f BED -p 0.01 --nomodel --shift 75 --extsize 150 -B --SPMR --keep-dup all --call-summits -n ${file::-4}
done

#We convert to a bed file as this is recommended by various people.
#Similar arguments to encode are used.
for file in *SPLIT*.bam
do
	bedtools bamtobed -i $file > ${file::-4}.bed
	macs2 callpeak -t ${file::-4}.bed -f BED -p 0.01 --nomodel --shift 75 --extsize 150 -B --SPMR --keep-dup all --call-summits -n ${file::-4}
done
