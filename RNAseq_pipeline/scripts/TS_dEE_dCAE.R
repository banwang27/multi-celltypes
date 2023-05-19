library(stringr)
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/EE/'
deseq2_path='/scratch/users/banwang3/celllines/data_backup/scratch/users/banwang3/celllines/analysis/DESeq2_norm/merged/'
celltypes=c('CM','HP','MN','PP','SKM')
countcutoff=25
dEE=c()
deseq2_all=c()
for (celltype in celltypes){
    df = read.csv(paste0(input_path,'RNA/dEE_zeroing_for',celltype,'_excludeRPE.txt'),sep='\t')
    deseq2 = read.csv(paste0(deseq2_path,celltype,'_merged_qvalue_DESeq2.txt'),sep='\t')
    idx = which(df$humallele<countcutoff & df$chpallele<countcutoff)
    df=df[-idx,]
    df$celltype=celltype  
    deseq2 = read.csv(paste0(deseq2_path,celltype,'_merged_qvalue_DESeq2.txt'),sep='\t') 
    deseq2=deseq2[,c('gene','log2FoldChange.hum','qval_mine.hum')]
    colnames(deseq2)=c('gene','RNA_lfc','RNA_FDR')
    deseq2$celltype=celltype
    dEE=rbind(dEE,df)
    deseq2_all=rbind(deseq2_all,deseq2)
}
dEE = merge(dEE,deseq2_all,by=c('gene','celltype'))
dCA=c()
for (celltype in celltypes){
    df = read.csv(paste0(input_path,'ATAC/dEE_for',celltype,'_includeSKM.txt'),sep='\t')
    df1 = read.csv(paste0('/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Filtered/',celltype,'_ATAC_Filtered_All_Peaks_fixed.txt'),sep='\t')
    df = df[,c('Peak.Humreffed','dEE.Humreffed')]
    colnames(df)=c('peak_id','dCA')
    df1 = df1[,c('Peak.Humreffed','l2fc_Humr','Human.binom_pval','Human.FDR')]
    colnames(df1)=c('peak_id','ATAC_lfc','ATAC_pval','ATAC_FDR')
    df2 = merge(df,df1)
    df2$celltype=celltype
    dCA=rbind(dCA,df2)
}
dCA = dCA[dCA$ATAC_pval<=0.05,]
peak_ids = unique(dCA$peak_id)
peak_ids=data.frame(peak_id = peak_ids,type=ifelse(grepl('promoter',peak_ids),'Promoter','Enhancer'))
dfopen_prom = peak_ids[peak_ids$type=='Promoter',]
dfopen_enha = peak_ids[peak_ids$type=='Enhancer',]
dfopen_prom$gene = str_split_fixed(dfopen_prom$peak_id,'promoter_',2)[,2]
while(length(which(grepl('_',dfopen_prom$gene)))>0){
    idx = which(grepl('_',dfopen_prom$gene))
     uniquegenedf = dfopen_prom[-idx,]
     dfopen_prom = dfopen_prom[idx,]
     gene1 = str_split_fixed(dfopen_prom$gene,'_',2)[,1] 
     gene2 = str_split_fixed(dfopen_prom$gene,'_',2)[,2] 
     tmp1 = dfopen_prom
     tmp1$gene = gene1
     tmp2 = dfopen_prom
     tmp2$gene = gene2
     dfopen_prom = rbind(uniquegenedf,tmp1,tmp2)
}
dfopen_enha$gene = str_split_fixed(dfopen_enha$peak_id,'enhancer_',2)[,2]
dfopen_enha$gene = str_split_fixed(dfopen_enha$gene,'_',2)[,1]
peak_ids = rbind(dfopen_prom,dfopen_enha)
dCA = merge(dCA,peak_ids)
final=merge(dCA,dEE,by=c('gene','celltype'))
final$samedirection = ifelse(final$ATAC_lfc*final$ASE>0,T,F)
dEEcutoff=0.75
dCAcutoff=0.75
idx = which(final$dEE>= dEEcutoff & final$dCA>= dCAcutoff)
final1 = final[idx,]
final1=final1[order(final1$celltype,-abs(final1$ATAC_lfc)),]
final2 = final1[final1$samedirection==T,]
write.table(final1,file=paste0('TS_list_ATAC_RNA_dEEgreaterorequalthan',dEEcutoff,'_dCAgreaterorequalthan',dCAcutoff,'.txt'),quote=F,sep='\t',row.names=F,col.names=T)
write.table(final2,file=paste0('TS_list_ATAC_RNA_dEEgreaterorequalthan',dEEcutoff,'_dCAgreaterorequalthan',dCAcutoff,'_samedirection.txt'),quote=F,sep='\t',row.names=F,col.names=T)
