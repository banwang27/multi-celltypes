#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G

#This script is essentially a way to do Hornet on a per chromosome basis.
samtools index rmdup.bam
samtools view -b rmdup.bam chr1 > rmdup.chr1.bam
samtools view -b rmdup.bam chr2 > rmdup.chr2.bam
samtools view -b rmdup.bam chr3 > rmdup.chr3.bam
samtools view -b rmdup.bam chr4 > rmdup.chr4.bam
samtools view -b rmdup.bam chr5 > rmdup.chr5.bam
samtools view -b rmdup.bam chr6 > rmdup.chr6.bam
samtools view -b rmdup.bam chr7 > rmdup.chr7.bam
samtools view -b rmdup.bam chr8 > rmdup.chr8.bam
samtools view -b rmdup.bam chr9 > rmdup.chr9.bam
samtools view -b rmdup.bam chr10 > rmdup.chr10.bam
samtools view -b rmdup.bam chr11 > rmdup.chr11.bam
samtools view -b rmdup.bam chr12 > rmdup.chr12.bam
samtools view -b rmdup.bam chr13 > rmdup.chr13.bam
samtools view -b rmdup.bam chr14 > rmdup.chr14.bam
samtools view -b rmdup.bam chr15 > rmdup.chr15.bam
samtools view -b rmdup.bam chr16 > rmdup.chr16.bam
samtools view -b rmdup.bam chr17 > rmdup.chr17.bam
samtools view -b rmdup.bam chr18 > rmdup.chr18.bam
samtools view -b rmdup.bam chr19 > rmdup.chr19.bam
samtools view -b rmdup.bam chr20 > rmdup.chr20.bam
samtools view -b rmdup.bam chr21 > rmdup.chr21.bam
samtools view -b rmdup.bam chr22 > rmdup.chr22.bam
samtools view -b rmdup.bam chrX > rmdup.chrX.bam
samtools view -b rmdup.bam chrY > rmdup.chrY.bam

python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr1.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr2.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr3.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr4.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr5.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr6.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr7.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr8.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr9.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr10.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr11.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr12.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr13.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr14.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr15.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr16.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr17.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr18.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr19.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr20.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr21.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chr22.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chrX.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/find_intersecting_snps_for_single.py -p rmdup.chrY.bam /scratch/users/astarr97/Myriad_Hybrid_ATAC/analysis/chimp/SNPs/wasp_split_species

bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr1.remap.fq1.gz -2 rmdup.chr1.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr1.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr2.remap.fq1.gz -2 rmdup.chr2.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr2.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr3.remap.fq1.gz -2 rmdup.chr3.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr3.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr4.remap.fq1.gz -2 rmdup.chr4.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr4.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr5.remap.fq1.gz -2 rmdup.chr5.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr5.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr6.remap.fq1.gz -2 rmdup.chr6.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr6.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr7.remap.fq1.gz -2 rmdup.chr7.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr7.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr8.remap.fq1.gz -2 rmdup.chr8.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr8.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr9.remap.fq1.gz -2 rmdup.chr9.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr9.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr10.remap.fq1.gz -2 rmdup.chr10.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr10.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr11.remap.fq1.gz -2 rmdup.chr11.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr11.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr12.remap.fq1.gz -2 rmdup.chr12.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr12.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr13.remap.fq1.gz -2 rmdup.chr13.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr13.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr14.remap.fq1.gz -2 rmdup.chr14.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr14.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr15.remap.fq1.gz -2 rmdup.chr15.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr15.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr16.remap.fq1.gz -2 rmdup.chr16.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr16.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr17.remap.fq1.gz -2 rmdup.chr17.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr17.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr18.remap.fq1.gz -2 rmdup.chr18.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr18.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr19.remap.fq1.gz -2 rmdup.chr19.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr19.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr20.remap.fq1.gz -2 rmdup.chr20.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr20.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr21.remap.fq1.gz -2 rmdup.chr21.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr21.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chr22.remap.fq1.gz -2 rmdup.chr22.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chr22.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chrX.remap.fq1.gz -2 rmdup.chrX.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chrX.out.bam
bowtie2 -X 2000 --very-sensitive-local -p 16 -x /scratch/users/astarr97/bowtie2_indices/bowtie2_gen_chimp.idx -1 rmdup.chrY.remap.fq1.gz -2 rmdup.chrY.remap.fq2.gz |  samtools view  -u - | samtools sort -n - > remapped_Aligned.chrY.out.bam

python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr1.to.remap.bam remapped_Aligned.chr1.out.bam rmdup.chr1.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr2A.to.remap.bam remapped_Aligned.chr2A.out.bam rmdup.chr2A.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr2B.to.remap.bam remapped_Aligned.chr2B.out.bam rmdup.chr2B.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr3.to.remap.bam remapped_Aligned.chr3.out.bam rmdup.chr3.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr4.to.remap.bam remapped_Aligned.chr4.out.bam rmdup.chr4.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr5.to.remap.bam remapped_Aligned.chr5.out.bam rmdup.chr5.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr6.to.remap.bam remapped_Aligned.chr6.out.bam rmdup.chr6.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr7.to.remap.bam remapped_Aligned.chr7.out.bam rmdup.chr7.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr8.to.remap.bam remapped_Aligned.chr8.out.bam rmdup.chr8.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr9.to.remap.bam remapped_Aligned.chr9.out.bam rmdup.chr9.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr10.to.remap.bam remapped_Aligned.chr10.out.bam rmdup.chr10.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr11.to.remap.bam remapped_Aligned.chr11.out.bam rmdup.chr11.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr12.to.remap.bam remapped_Aligned.chr12.out.bam rmdup.chr12.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr13.to.remap.bam remapped_Aligned.chr13.out.bam rmdup.chr13.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr14.to.remap.bam remapped_Aligned.chr14.out.bam rmdup.chr14.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr15.to.remap.bam remapped_Aligned.chr15.out.bam rmdup.chr15.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr16.to.remap.bam remapped_Aligned.chr16.out.bam rmdup.chr16.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr17.to.remap.bam remapped_Aligned.chr17.out.bam rmdup.chr17.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr18.to.remap.bam remapped_Aligned.chr18.out.bam rmdup.chr18.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr19.to.remap.bam remapped_Aligned.chr19.out.bam rmdup.chr19.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr20.to.remap.bam remapped_Aligned.chr20.out.bam rmdup.chr20.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr21.to.remap.bam remapped_Aligned.chr21.out.bam rmdup.chr21.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chr22.to.remap.bam remapped_Aligned.chr22.out.bam rmdup.chr22.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chrX.to.remap.bam remapped_Aligned.chrX.out.bam rmdup.chrX.kept.bam
python /scratch/users/astarr97/Myriad_Hybrid_ATAC/Hornet/mapping/filter_remapped_reads.py -p rmdup.chrY.to.remap.bam remapped_Aligned.chrY.out.bam rmdup.chrY.kept.bam

samtools sort -o rmdup.chr1.kept.sort.bam rmdup.chr1.kept.bam
samtools sort -o rmdup.chr2A.kept.sort.bam rmdup.chr2A.kept.bam
samtools sort -o rmdup.chr2B.kept.sort.bam rmdup.chr2B.kept.bam
samtools sort -o rmdup.chr3.kept.sort.bam rmdup.chr3.kept.bam
samtools sort -o rmdup.chr4.kept.sort.bam rmdup.chr4.kept.bam
samtools sort -o rmdup.chr5.kept.sort.bam rmdup.chr5.kept.bam
samtools sort -o rmdup.chr6.kept.sort.bam rmdup.chr6.kept.bam
samtools sort -o rmdup.chr7.kept.sort.bam rmdup.chr7.kept.bam
samtools sort -o rmdup.chr8.kept.sort.bam rmdup.chr8.kept.bam
samtools sort -o rmdup.chr9.kept.sort.bam rmdup.chr9.kept.bam
samtools sort -o rmdup.chr10.kept.sort.bam rmdup.chr10.kept.bam
samtools sort -o rmdup.chr11.kept.sort.bam rmdup.chr11.kept.bam
samtools sort -o rmdup.chr12.kept.sort.bam rmdup.chr12.kept.bam
samtools sort -o rmdup.chr13.kept.sort.bam rmdup.chr13.kept.bam
samtools sort -o rmdup.chr14.kept.sort.bam rmdup.chr14.kept.bam
samtools sort -o rmdup.chr15.kept.sort.bam rmdup.chr15.kept.bam
samtools sort -o rmdup.chr16.kept.sort.bam rmdup.chr16.kept.bam
samtools sort -o rmdup.chr17.kept.sort.bam rmdup.chr17.kept.bam
samtools sort -o rmdup.chr18.kept.sort.bam rmdup.chr18.kept.bam
samtools sort -o rmdup.chr19.kept.sort.bam rmdup.chr19.kept.bam
samtools sort -o rmdup.chr20.kept.sort.bam rmdup.chr20.kept.bam
samtools sort -o rmdup.chr21.kept.sort.bam rmdup.chr21.kept.bam
samtools sort -o rmdup.chr22.kept.sort.bam rmdup.chr22.kept.bam
samtools sort -o rmdup.chrX.kept.sort.bam rmdup.chrX.kept.bam
samtools sort -o rmdup.chrY.kept.sort.bam rmdup.chrY.kept.bam

samtools index rmdup.chr1.kept.sort.bam
samtools index rmdup.chr2A.kept.sort.bam
samtools index rmdup.chr2B.kept.sort.bam
samtools index rmdup.chr3.kept.sort.bam
samtools index rmdup.chr4.kept.sort.bam
samtools index rmdup.chr5.kept.sort.bam
samtools index rmdup.chr6.kept.sort.bam
samtools index rmdup.chr7.kept.sort.bam
samtools index rmdup.chr8.kept.sort.bam
samtools index rmdup.chr9.kept.sort.bam
samtools index rmdup.chr10.kept.sort.bam
samtools index rmdup.chr11.kept.sort.bam
samtools index rmdup.chr12.kept.sort.bam
samtools index rmdup.chr13.kept.sort.bam
samtools index rmdup.chr14.kept.sort.bam
samtools index rmdup.chr15.kept.sort.bam
samtools index rmdup.chr16.kept.sort.bam
samtools index rmdup.chr17.kept.sort.bam
samtools index rmdup.chr18.kept.sort.bam
samtools index rmdup.chr19.kept.sort.bam
samtools index rmdup.chr20.kept.sort.bam
samtools index rmdup.chr21.kept.sort.bam
samtools index rmdup.chr22.kept.sort.bam
samtools index rmdup.chrX.kept.sort.bam
samtools index rmdup.chrY.kept.sort.bam

samtools merge rmdup.remap.kept.bam rmdup.chr1.kept.sort.bam rmdup.chr2A.kept.sort.bam rmdup.chr2B.kept.sort.bam rmdup.chr3.kept.sort.bam rmdup.chr4.kept.sort.bam rmdup.chr5.kept.sort.bam rmdup.chr6.kept.sort.bam rmdup.chr7.kept.sort.bam rmdup.chr1.kept.sort.bam rmdup.chr8.kept.sort.bam rmdup.chr1.kept.sort.bam rmdup.chr9.kept.sort.bam rmdup.chr10.kept.sort.bam rmdup.chr11.kept.sort.bam rmdup.chr12.kept.sort.bam rmdup.chr13.kept.sort.bam rmdup.chr14.kept.sort.bam rmdup.chr15.kept.sort.bam rmdup.chr16.kept.sort.bam rmdup.chr17.kept.sort.bam rmdup.chr18.kept.sort.bam rmdup.chr19.kept.sort.bam rmdup.chr20.kept.sort.bam rmdup.chr21.kept.sort.bam rmdup.chr22.kept.sort.bam rmdup.chrX.kept.sort.bam rmdup.chrY.kept.sort.bam