# Description of files

ATAC_All_Peaks_Down_Gnomad_Constraint_Summary_Filtered.txt: 
gnomAD constraint metrics for each downsampled peak after removing all protein coding exons from the human gencode gtf file.

ATAC_All_Peaks_Gnomad_Constraint_Summary_Filtered.txt: 
the same except that it is for non-downsampled peaks.

Correct_Gene_Names.txt: 
map the Horlbeck et al. gene names used for peaks to the human gencode v38 gene names.

Differential expression/: 
files for differential expression analysis on publicly available bulk RNAseq data.

Differential_Expression_For_Blake_Pavlovic.R: 
differential expression analysis with DESeq2 for the Blake et al. and Pavlovic et al. data.It uses files like "Blake_Gilad_2020_Config_HumChp_Heart_nocovs.txt" to configure the script and the gene_length files here to compute TPM.

Ma_Etal_DESeq2_Config_Files/: 
config files for the pseudobulked Ma et al. DLPFC data.

Filtering_ATAC_And_Striatal_Organoid-Analysis.ipynb: 
filter and aggregate ATAC counts and analyse the striatal organoid data (mostly just compute CPM from it).

Kanton_Etal_Analysis.ipynb: 
pseudotime analysis applied to the Kanton et al. single cell data from the ventral forebrain interneuron lineage.

Processing_Ma_Etal_DLPFC_2022.ipynb: 
psuedobulking the DLPFC data from Ma et al.

Selection_Test_Z_Score.ipynb: 
perform tests for lineage-specific selection. Additional files needed to run this script can be found in the Zotero repo associated with Starr et al. 2023. 
