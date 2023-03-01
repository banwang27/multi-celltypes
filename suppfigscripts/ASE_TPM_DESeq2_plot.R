library(ggplot2)  
library(patchwork)
library(stringr)
library(dplyr)
# Which genes you want to plot?
target_genes=c('NTF3','NTRK1','NTRK3','SCN1B','SCN2B','SYT2','NR1H3')

#TPM plot with hyb label
pboxplot<-function(sampleTPM,gene){
          p<-ggplot(data=sampleTPM,aes(x=celltype,y=TPM))+
          geom_boxplot(outlier.size=0,fill="white",colour="black",width=0.5)+
        geom_point(data=sampleTPM,aes(x=celltype,y=TPM,fill=hyb), size=1.5, shape=21, colour="grey20",position=position_jitter(width=0.1, height=0.1))+
          xlab('')+  
          ylab(paste0(gene, ' TPM'))+
    theme_bw()+  
          guides(fill=guide_legend(title="Hybrid line"))+
          theme(legend.position="top",panel.grid.major=element_blank(),panel.grid.minor=element_blank())
      return(p)
}
#Allele plot with Allele label
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
#DEgene plot with siginificance level label
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
celltypes=c('CM','HP','MN','PP','RPE','SKM')
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
df=read.csv(paste0(backup_path,input_path,'TPM.txt'),sep='\t')
idx =grep('total',colnames(df))
idx1=grep('gene',colnames(df))
df = df[,c(idx1,idx)]
tmp=c()
for (celltype in celltypes){
     idx = grep(celltype,colnames(df))
     tmp = cbind(tmp,rowMeans(df[,idx]))
}
colnames(tmp)=celltypes
tmp = data.frame(tmp)
tmp$gene=df$gene
df1 = df[df$gene %in% target_genes,]
genenames=df1$gene
idx = grep('gene',colnames(df1))
sampleTPM = c()
for (gene in genenames){
     tmp1 = t(df1[df1$gene==gene,-idx]) 
     tmp1=as.data.frame(tmp1)
     colnames(tmp1)='TPM' 
      tmp1$sampleid=row.names(tmp1)
      tmp1$gene=gene
      sampleTPM=rbind(sampleTPM,tmp1)
}
sampleTPM$celltype=str_split_fixed(sampleTPM$sampleid,'_',4)[,2]
sampleTPM$hyb=str_split_fixed(sampleTPM$sampleid,'_',4)[,3]
sampleTPM$celltype=factor(sampleTPM$celltype)
result=deseq2all[deseq2all$gene %in% target_genes,]
result$plab<-ifelse(result$qval_mine<0.001,'***',ifelse(result$qval_mine<0.01,'**',ifelse(result$qval_mine<0.05,'*','')))
for (gene in target_genes){
    TPMplot<-pboxplot(sampleTPM[sampleTPM$gene==gene,],gene)
    alleleplot<-palleleplot(final[final$gene==gene,],gene)
    deseq2plot<-pdeseq2plot(result[result$gene==gene,],gene)
    p<-TPMplot|alleleplot|deseq2plot
    pdf(paste0(output_path,'TPM_Allele_Deseq2_',gene ,'.pdf'),width=12,height=3.5, useDingbats=FALSE)
    print(p)
    dev.off()
}
