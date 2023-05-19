library(ggplot2)  
library(stringr)
library(dplyr)
library(patchwork)
align2='human'
TPM = read.csv('/scratch/users/banwang3/celllines/analysis/ASE_norm/human/filtered/rawcounts.txt',sep='\t')
idx =grep('total',colnames(TPM))
TPM = TPM[,-idx]
idx = grep('RNA',colnames(TPM))
CPM = sweep(as.matrix(TPM[,idx]), 2, as.double(colSums(TPM[,idx])/1000000), `/`)
CPM = data.frame(CPM)
CPM$gene = TPM$gene
gene = 'CTSF'
df = data.frame(t(CPM[CPM$gene==gene,-ncol(CPM)]))
colnames(df)=c('CPM')    
df$sample_name = row.names(df)
df$celltype = str_split_fixed(df$sample_name,'_',4)[,2]
df$Hyb = str_split_fixed(df$sample_name,'_',4)[,3]
df$Allele = str_to_title(str_split_fixed(df$sample_name,'_',4)[,4])
df$sample = str_split_fixed(df$sample_name,'_',4)[,1]
df$Allele = ifelse(df$Allele == 'Hum','Human','Chimp')
df$Allele = factor(df$Allele, levels=c('Human','Chimp'))
df$celltype = factor(df$celltype)
df = df[df$celltype!='RPE',]
df1 =df%>% group_by(celltype,Allele)%>%summarise(CPM=mean(CPM)) 
df1 = as.data.frame(df1)
df1$Allele = factor(df1$Allele, levels=c('Human','Chimp'))
df1$celltype = factor(df1$celltype)
barwidth = 0.6
pointsize=0.3

celltypes = c('CM','HP','MN','PP','SKM')
alleleplot<-function(df,celltype){
    tmpdf = df[df$celltype==celltype,]
    tmpdf1 = as.data.frame(tmpdf%>%group_by(Allele)%>%summarise(CPM=mean(CPM)))
    tmpdf1$Allele =factor(tmpdf1$Allele, levels=c('Chimp','Human'))
    p1<-ggplot()+
        geom_bar(data=tmpdf1,aes(x=Allele,y=CPM,fill=Allele),colour="black",stat="identity",width=barwidth, alpha=0.8)+ 
        theme_bw()+ 
        xlab('')+
        ylab('')+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        scale_fill_manual(values=c("#1E90FF","#FF6A6A"))+
        coord_flip()+
        theme(axis.text.x = element_text(size = 8),legend.position="none") 
    return(p1)
}
for (celltype in celltypes){
    p<-alleleplot(df,celltype)
    pdf(paste0('/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/Allele_CPM_',gene,'_',celltype,'.pdf'),width=2,height=1.5, useDingbats=FALSE)
    print(p)
    dev.off()
}

