library(ggplot2)  
library(patchwork)
library(stringr)
library(dplyr)
palleleplot<-function(sampleCPM,gene){
      p<-ggplot(data=sampleCPM,aes(x=celltype,y=CPM,fill=Allele))+
          geom_bar(position = "dodge",colour="black",stat="summary",fun="mean",width=0.6,  alpha = 0.8) +                                                                   
          geom_point(data=sampleCPM,aes(x=celltype,y=CPM,fill=Allele), size=1.5, shape=21, colour="grey20",position = position_dodge(width=0.6),show.legend=F)+
          ylab(paste0(gene,' CPM'))+
          xlab('')+  
          theme_bw()+                                      
          scale_fill_manual(values=c("#FF6A6A","#1E90FF"))+                                                                                      
          theme(axis.title = element_text(size = 8),legend.position="top",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 return(p)                                                        
}       
pdeseq2plot<-function(result,gene){
    p<-ggplot(data=result,aes(x=celltype,y=log2FoldChange))+
        geom_hline(yintercept=0)+
        geom_bar(stat="identity",width=0.5,alpha=0.8)+
        geom_text(aes(x=celltype,label=plab),vjust="inward")+
        ylab(paste0(gene,' Log2(Human/Chimp)'))+
        xlab('')+
         theme_bw()+ 
         theme(axis.title = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    return(p)
}
celltypes=c('CM','HP','MN','PP','SKM')
backup_path='/scratch/users/banwang3/celllines/data_backup'
input_path='/scratch/users/banwang3/celllines/analysis/ASE_norm/human/filtered/'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/suppplots/' 
deseq2_path='/scratch/users/banwang3/celllines/analysis/DESeq2_norm/human/'
endfile='_aligned2human_qvalue_DESeq2.txt'
deseq2all=c()
for (celltype in celltypes){
    deseq2=read.csv(paste0(backup_path,deseq2_path,celltype,endfile),sep='\t')
    deseq2=deseq2[,c('gene','log2FoldChange','qval_mine')]
    deseq2$celltype=celltype
    deseq2all=rbind(deseq2all,deseq2)
}
TPM = read.csv(paste0(backup_path,input_path,'rawcounts.txt'),sep='\t')
idx =grep('total',colnames(TPM))
TPM = TPM[,-idx]
idx = grep('RNA',colnames(TPM))
CPM = sweep(as.matrix(TPM[,idx]), 2, as.double(colSums(TPM[,idx])/1000000), `/`)
CPM = data.frame(CPM)
CPM$gene = TPM$gene
target_genes=c('GAD1')
#target_genes=c('FABP7')
gene=target_genes[1]
final=c()
for (gene in target_genes){
    df = data.frame(t(CPM[CPM$gene==gene,-ncol(CPM)]))
    colnames(df)=c('CPM')    
    df$sample_name = row.names(df)
    df$celltype = str_split_fixed(df$sample_name,'_',4)[,2]
    df$Hyb = str_split_fixed(df$sample_name,'_',4)[,3]
    df$Allele = str_split_fixed(df$sample_name,'_',4)[,4]
    df$sample = str_split_fixed(df$sample_name,'_',4)[,1]
    df$Allele = ifelse(df$Allele == 'hum','Human','Chimp')
    df$gene=gene
    final=rbind(final,df)
}
final$Allele=factor(final$Allele, levels=c('Human','Chimp'))
final$celltype=factor(final$celltype)
final=final[final$celltype!='RPE',]
meanCPM = data.frame(final%>% group_by(gene,celltype,Allele)%>% summarise(CPM=mean(CPM)))     
result=deseq2all[deseq2all$gene %in% target_genes,]
result$plab<-ifelse(result$qval_mine<0.001,'***',ifelse(result$qval_mine<0.01,'**',ifelse(result$qval_mine<0.05,'*','')))

alleleplot<-palleleplot(final[final$gene==gene,],gene)
deseq2plot<-pdeseq2plot(result[result$gene==gene,],gene)
p<-alleleplot/deseq2plot
  
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Unfiltered/'
peak_id='peak_CM_HP_MN_SKM_187930_enhancer_GAD1_6614.5_ERICH2_39630.5'
#peak_id='peak_HP_MN_PP_SKM_295546_promoter_FABP7'
ATAC=c()
lfc=c()
for (celltype in celltypes){
    df = read.csv(paste0(input_path,celltype,'_ATAC_Filtered_All_Peaks_fixed_nofilt.txt'),sep='\t')    
    idx = grep('humreffed',colnames(df))
    df1 = as.data.frame(sweep(as.matrix(df[,idx]), 2, as.double(colSums(df[,idx])/1000000), `/`))
    lfc = c(lfc,df[df$Peak.Humreffed==peak_id,'l2fc_Humr'])
    df1$peak_id=df$Peak.Humreffed
    df2 = df1[df1$peak_id==peak_id,-ncol(df1)]
    df2 = data.frame(t(df2))
    colnames(df2)='CPM' 
    df2$Allele=ifelse(grepl('Human',row.names(df2)),'Human','Chimp')
    df2$Hyb=ifelse(grepl('Hyb1',row.names(df2)),'Hyb1','Hyb2')
    df2$celltype=celltype
    ATAC = rbind(ATAC,df2)
}
ATAC$Allele=factor(ATAC$Allele,levels=c('Human','Chimp'))
lfc = data.frame(lfc=lfc,celltype=celltypes)
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Filtered/'
pval=c()
FDR=c()
l2fc=c()
for (celltype in celltypes){
    df = read.csv(paste0(input_path,celltype,'_ATAC_Filtered_All_Peaks_fixed.txt'),sep='\t')
    idx = which(df$Peak.Humreffed==peak_id)
    if (length(idx)>0){
        pval=c(pval,df[idx,'Human.binom_pval'])
        FDR=c(FDR,df[idx,'Human.FDR'])
        l2fc=c(l2fc,df[idx,'l2fc_Humr'])
    }else{
        pval=c(pval,1)
        FDR=c(FDR,1)
        l2fc=c(l2fc,0)
    }
}
data = data.frame(celltype=celltypes,pval=pval,FDR=FDR,l2fc=l2fc)
lfc = merge(lfc,data)
lfc$plab=ifelse(lfc$pval<0.001,'***',ifelse(lfc$pval<0.01,'**',ifelse(lfc$pval<0.05,'*','')))
colnames(lfc)=c('celltype','log2FoldChange','pval','FDR','l2fc','plab')
alleleplot<-palleleplot(ATAC,gene)
deseq2plot<-pdeseq2plot(lfc,gene)
p1<-alleleplot/deseq2plot
p2<-p1|p
pdf(paste0(output_path,'RNA_ATAC_Allele_lfc_',gene ,'.pdf'),width=8,height=6.5, useDingbats=FALSE)
print(p2)
dev.off()
