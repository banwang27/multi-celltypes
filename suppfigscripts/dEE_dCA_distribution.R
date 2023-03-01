library(ggplot2)
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/EE/'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/suppplots/'
keyword='RNA'
includeRPE=T
if (includeRPE){
    celltypes=c('CM','HP','MN','PP','RPE','SKM')
    endname='includeSKM'
}else{
    celltypes=c('CM','HP','MN','PP','SKM')
    endname='excludeRPE'
}
countcutoff=25
final=c()
for (celltype in celltypes){
    df = read.csv(paste0(input_path,keyword,'/dEE_zeroing_for',celltype,'_',endname,'.txt'),sep='\t')
    idx = which(df$humallele<countcutoff & df$chpallele<countcutoff)
    df=df[-idx,]
    df$celltype=celltype
    final=rbind(final,df)
}
final$ASElab=ifelse(abs(final$ASE)>=0.5,'absLFC >= 0.5','absLFC < 0.5')
p<-ggplot(final, aes(x=dEE,fill = ASElab)) + 
    geom_histogram(bins=100)+
    facet_wrap(~celltype,scales="free_y")+
    theme_bw()+
    guides(fill=guide_legend(title=""))+
    theme(legend.position="top",panel.grid.major=element_blank(),panel.grid.minor=element_blank())
pdf(paste0(output_path,'dEE_',keyword,'_distribution_',endname,'.pdf'),useDingbats=FALSE,width=6,height=5) 
print(p)
dev.off()
keyword='ATAC'
celltypes=c('CM','HP','MN','PP','SKM')
final=c()
SKM='includeSKM'
for (celltype in celltypes){
    df = read.csv(paste0(input_path,keyword,'/dEE_for',celltype,'_',SKM,'.txt'),sep='\t')
    df1 = read.csv(paste0('/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Filtered/',celltype,'_ATAC_Filtered_All_Peaks_fixed.txt'),sep='\t')
    df = df[,c('Peak.Humreffed','dEE.Humreffed')]
    colnames(df)=c('peak_id','dEE')
    df1 = df1[,c('Peak.Humreffed','l2fc_Humr','Human.binom_pval','Human.FDR')]
    colnames(df1)=c('peak_id','lfc','pval','FDR')
    df2 = merge(df,df1)
    df2$celltype=celltype
    final=rbind(final,df2)
}
final$ASElab=ifelse(abs(final$lfc)>=0.5,'absLFC >= 0.5','absLFC < 0.5')
p<-ggplot(final, aes(x=dEE,fill = ASElab)) + 
    geom_histogram(bins=100)+
    facet_wrap(~celltype,scales="free_y")+
    theme_bw()+
    guides(fill=guide_legend(title=""))+
    theme(legend.position="top",panel.grid.major=element_blank(),panel.grid.minor=element_blank())
pdf(paste0(output_path,'dEE_',keyword,'_distribution.pdf'),useDingbats=FALSE,width=6,height=5)
print(p) 
dev.off()

