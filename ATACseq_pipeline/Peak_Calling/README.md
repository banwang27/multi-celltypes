# Steps for peak calling
Read `pipeline.sh` for more details.

The general pipeline we follow is this:
1. Map reads to the genome with `bowtie2`, then correct for mapping bias with `Hornet` (Bowtie, Hornet and SAM tools, see ../Mapping).
2. Filter out any multimapping reads (`Picard`).
3. Filter out any reads mapping outside of chrX and autosomes (`samtools`).
4. Use `macs2` to call peaks on a merged bam file containing all the reads for a particular condition as well as the individual bam files (`macs2`).
5. Filter out any peaks that were not called in all replicates that went into a merged bam file (awk, in `pipeline.sh`).
6. Perform liftover on the narrowpeak files such that we lift human to chimp and back to human and lift chimp to human (`liftOver`).
7. Intersect both now human referenced files (`bedtools`).
8. Go through the intersected file and merge the peaks such that we take the union of any overlapping intervals whether they come from the chimpanzee referenced file or the human referenced file (`find_partners.py`, `merge.py`).
This outputs a narrowPeak file in the traditional format as well as a partners file that identifies homologous peaks.
9. Intersect with the Encode blacklist regions and filter those on out (`bedtools` and `awk`, see `pipeline.sh`).
10. Merge any peaks with overlap (`merge.py`).
11. Lift the human referenced file back to the chimpanzee reference and merge any overlapping peaks there (`liftOver`, `merge.py`).
12. Filter out any peaks that, as a result of merging, no longer map one to one between the human and chimpanzee reference genomes using the partners file (`filter_11.py`).
13. Merge with promoters (again taking the union) and annotate any non-promoters as enhancers attached to the two closest genes (`bedtools`, `anno_prom.py`, `anno_enh.py`).
14. Take this and run through the modified pipeline in ../Merging_Peaks to create a unified peaklist across all cell types.
15. Convert to a gtf file (`np_to_gtf.py`).
16. Split the bam files into human and chimpanzee reads (see ../Counting).
17. Count reads in these bam files with HTSeq (see ../Counting).

Processing promoters:
1. Downloaded the list of TSS from https://elifesciences.org/articles/19760.  
2. Expanded +/- 1000 bp around the TSS and merge (with the union) any overlapping promoters and lifted over to hg38.
3. Lifted them over to chimp and back to human and filtered out any that failed to liftover from the original list as well as the chimpanzee list.

The list of promoters is included here:

`Human_Promoters_Ortho_Sorted_hg38_PanTro6.bed`;
`Human_Promoters_Ortho_Sorted_PanTro6_PanTro6.bed`

Note:
If there is only one replicate, we need to split the bamfile. 
To accomplish this, we use `prep_1_rep.sh`.
This uses Picard to sort the bamfile (must use Picard to sort the bamfile) and then Picard again to split the bamfile.  
After that, we treat the two split bamfiles as pseudoreplicates and go through the above pipeline.
