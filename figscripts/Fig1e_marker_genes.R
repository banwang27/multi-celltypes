library(ggplot2)
library(stringr)
backup_path='/scratch/users/banwang3/celllines/data_backup'
df = read.csv(paste0(backup_path,'/scratch/users/banwang3/celllines/analysis/ASE_norm/human/filtered/TPM.txt'),sep='\t')
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/'
idx =grep('total',colnames(df))
idx1=grep('gene',colnames(df))
df = df[,c(idx1,idx)]
celltypes = c('CM','HP','PP','SKM','MN','RPE')
tmp=c()
for (celltype in celltypes){
    idx = grep(celltype,colnames(df))
    tmp = cbind(tmp,rowMeans(df[,idx]))
}
colnames(tmp)=celltypes
tmp = data.frame(tmp)
tmp$gene=df$gene
target_genes=c('NKX2-5','HNF1B','OLIG2','GP2','PAX6','MYLPF')
df1 = df[df$gene %in% target_genes,]
genenames=df1$gene
target_genes=genenames
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
sampleTPM$hyb=ifelse(sampleTPM$hyb=='Hyb1','Hybrid1','Hybrid2')
pboxplot<-function(sampleTPM){
    p<-ggplot(data=sampleTPM,aes(x=celltype,y=TPM))+
        geom_boxplot(outlier.size=0,fill="white",colour="black",width=0.5)+
        geom_point(data=sampleTPM,aes(x=celltype,y=TPM,fill=hyb), size=1.5, shape=21, colour="grey20",position=position_jitter(width=0.1, height=0.1))+
        scale_fill_manual(values=c("#6959cd","#3cb371")) +
        facet_wrap(~gene,scales="free")+
        xlab('')+
        ylab('Expression level (TPM)')+
        theme_bw()+    
        guides(fill=guide_legend(title="Hybrid line"))+
        theme(legend.position="top",panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    return(p)
}
final_genes=c('NKX2-5','HNF1B','OLIG2','GP2','PAX6','MYLPF')
sampleTPM1 = sampleTPM[sampleTPM$gene %in% final_genes,]
sampleTPM1$gene=factor(sampleTPM1$gene, levels=final_genes)
pdf(paste0(output_path,'marker_gene_TPM.pdf'),useDingbats=FALSE,width=7,height=5)
p<-pboxplot(sampleTPM1)
print(p)
dev.off()
