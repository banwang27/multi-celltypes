peak_id="peak_HP_MN_PP_SKM_295546_promoter_FABP7"
peak_short="peak_HP_MN_PP_SKM_295546_"
folder_name="FABP7"


ml python/3.6.1
ml R/4.2.0
cd /scratch/users/banwang3/references/outputs/
mkdir ${folder_name}
cd /scratch/users/banwang3/references/outputs/${folder_name}
# Get Peak coordinate
grep ${peak_short} /scratch/users/banwang3/celllines/ATAC/Peaklist_Final/All_Peaks_Peaklist_Final_Anno_Humreffed_exchr20.bed > hg38_peak_coord.bed
grep ${peak_short} /scratch/users/banwang3/celllines/ATAC/Peaklist_Final/All_Peaks_Peaklist_Final_Anno_Chpreffed_exchr20.bed > pantro6_peak_coord.bed

#Extend to 4096bp window
Rscript /scratch/users/banwang3/references/scripts/extendto4096bpwindow.R --generator ${folder_name}

bedtools  getfasta -fi /scratch/users/banwang3/references/human/human.fasta -bed hg38_centeratpeak_4096bp.bed -fo hg38_centeratpeak_4096bp.fasta
bedtools  getfasta -fi /scratch/users/banwang3/references/chimp/chimp.fasta -bed pantro6_centeratpeak_4096bp.bed -fo pantro6_centeratpeak_4096bp.fasta

bedtools intersect -a /scratch/users/banwang3/references/human/ASE_SNPs.FILTER.SPLIT_SPECIES_HUMAN.bed -b hg38_centeratpeak_4096bp.bed -wa -wb > SNPlist_hg38_4096bp.bed
bedtools intersect -a /scratch/users/banwang3/references/chimp/ASE_SNPs.FILTER.SPLIT_SPECIES_CHIMP.bed -b pantro6_centeratpeak_4096bp.bed -wa -wb > SNPlist_pantro6_4096bp.bed


CrossMap.py bed /scratch/users/banwang3/celllines/liftover/hg38ToPanTro6.over.chain.gz SNPlist_hg38_4096bp.bed SNPlist_hg38_4096bp_lift2pantro6.bed
CrossMap.py bed /scratch/users/banwang3/celllines/liftover/panTro6ToHg38.over.chain.gz SNPlist_pantro6_4096bp.bed SNPlist_pantro6_4096bp_lift2hg38.bed

Rscript /scratch/users/banwang3/references/scripts/make4096bpbed_centeratSNP.R --generator ${folder_name}

filenames="hg38_centeratpeak_4096bp.bed hg38_centeratSNP_4096bp.bed pantro6_centeratSNP_4096bp_lift2hg38.bed"
for filename in $filenames;
do
    bedtools  getfasta -fi /scratch/users/banwang3/references/human/human.fasta -bed ${filename} -fo ${filename::-4}.fasta
done
filenames2="pantro6_centeratpeak_4096bp.bed pantro6_centeratSNP_4096bp.bed hg38_centeratSNP_4096bp_lift2pantro6.bed"
for filename in $filenames2;
do
    bedtools  getfasta -fi /scratch/users/banwang3/references/chimp/chimp.fasta -bed ${filename} -fo ${filename::-4}.fasta
done

Rscript /scratch/users/banwang3/references/scripts/base_changer.R --generator ${folder_name}

cat hg38_centeratSNP_4096bp.fasta hg38_centeratSNP_4096bp_lift2pantro6.fasta SNPfasta_hg38_centeratSNP_4096bp.fasta > hg38_all.fasta
cat pantro6_centeratSNP_4096bp_lift2hg38.fasta pantro6_centeratSNP_4096bp.fasta SNPfasta_pantro6_centeratSNP_4096bp.fasta > pantro6_all.fasta
