#!/bin/bash
#SBATCH --job-name=bamcoverage
#SBATCH --time=3-00:00:00
##SBATCH --qos=long
#SBATCH --mem=16GB
#SBATCH -p hbfraser
#SBATCH --mail-type=END,FAIL


ml python/3.6.1
ml biology
ml samtools
ml py-deeptools
cd /scratch/users/banwang3/celllines/ATAC/human/bam/
celltypes="CM HP MN PP SKM"
# bamCoverage bam to bigwig and normalizing on binSize=1
for bamfile in *humreffed.merged.Chpe.sort.bam;
do
    samtools view -L /scratch/users/banwang3/references/human/exchr20.bed -o bigwig_1bp/${bamfile::-4}.bam $bamfile
    samtools index bigwig/${bamfile::-4}.bam
    echo transfering $bamfile
    bamCoverage -b $bamfile -o bigwig_1bp/${bamfile::-4}.bw --binSize 1 --normalizeUsing CPM --effectiveGenomeSize 2862010578 --ignoreForNormalization chr20 --extendReads
done

# Use bigwigCompare to compare two tracks
for celltype in $celltypes;
do
    echo transfering $celltype
    bigwigCompare -b1 Hume/bigwig_1bp/ATAC_${celltype}_humreffed.merged.Hume.sort.bw -b2 Chpe/bigwig_1bp/ATAC_${celltype}_humreffed.merged.Chpe.sort.bw --pseudocount 1 --skipZeroOverZero --operation log2 -bs 1 -r chr11:4696850:4698850 -of bedgraph -o ${celltype}_Hume_Chpe_bigwigCompare_OR51E2.bedgraph
done


# Use pyGenomeTracks to visualize tracks
pyGenomeTracks --tracks Fig3c_CTSF_both_allele.ini --region chr11:66540000-66590000 -o CTSF_tracks.pdf
