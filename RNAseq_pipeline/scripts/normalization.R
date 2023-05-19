library(stringr)
align2='human'
input_path=paste0('/scratch/users/banwang3/celllines/analysis/ASE/',align2,'/')
output_path=paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',align2,'/')
humangenes=read.csv('/scratch/users/banwang3/celllines/scripts/RScripts/All_Human_Genes.csv',sep=',')
samplenames<-c('RNA3_CM_Hyb1','RNA4_CM_Hyb2')
endfile='_ase_by_reads_merged.txt'
# Filter out reads on chr20 and pseudogenes, then perform normalization to balance human and chip allele ratios
for (samplename in samplenames){
    df = read.csv(paste0(input_path,samplename,endfile),sep='\t',skip=1)
    df1<-df[df$chrom!='chr20',]
    refalleles<-sum(df1$ref_counts)
    altalleles<-sum(df1$alt_counts)
    ratio<-refalleles/altalleles
    if (ratio<1){   
        new_alt = round(df1$alt_counts*ratio)
        df1$alt_counts<-new_alt
    }else{
        new_ref = round(df1$ref_counts*(1/ratio))
        df1$ref_counts<-new_ref
    }
    filename = paste0(output_path,samplename,'_ase_normalized.txt')
    write.table(df1,filename,sep='\t',quote=F,row.names=F,col.names=T)
    df2 = df1[df1$gene %in% humangenes$Gene_Symbol,]         
    write.table(df2,paste0(output_path,'filtered/',samplename,'_ase_normalized.txt'),sep='\t',quote=F,row.names=F,col.names=T) 
}
# Filter out reads on chr20 and pseudogenes
for (samplename in samplenames){
    df = read.csv(paste0(input_path,samplename,endfile),sep='\t',skip=1)
    df1<-df[df$chrom!='chr20',]
    df2<-df1[df1$gene %in% humangenes$Gene_Symbol,]
    write.table(df2,paste0(input_path,'filtered/',samplename,endfile),sep='\t',quote=F,row.names=F,col.names=T)
}
