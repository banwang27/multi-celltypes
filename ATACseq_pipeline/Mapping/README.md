# Steps for mapping 

To map ATAC reads, first set up the directory with the Hornet and Scripts folders, the SNPs in the analysis folder, picard.jar, the attribute file etc.

Next, use Snakefile and runSnake.sh to get to the rmdup.bam stage (bam file with duplicates removed).

Next, use human_by_chr.sh and chimp_by_chr.sh in each folder to do the hornet process.

after_by_chr.sh then finishes this and produces the final file for peak calling.

Please contact Fraser Lab at Stanford University for resulting SNP lists.

This process is described in more detail in Agoglia et al. 2021.
