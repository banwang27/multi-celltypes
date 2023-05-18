library(DESeq2)
library(stringr)
library(ggplot2)
Aligned2='human'
all_files = list.files(paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',Aligned2,'/'))
sub_files = all_files[grep('_ase_normalized.txt', all_files)]
#sub_files=sub_files[-grep(cell_type,sub_files)]
input_path <- paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',Aligned2,'/')
endFile <- '_ase_normalized.txt'
output_path <-paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',Aligned2,'/deseq2_normalized_counts/')
valMat = c()
#allhumangenes<-read.table('/scratch/users/banwang3/celllines/analysis/ASE/human/withPseudogenes/All_Human_Genes.csv',sep=',',header=TRUE)
for (a in sub_files){  
    input_filename<-paste0(input_path,a)
    repVals <- read.table(input_filename , header=TRUE, na.strings = "", sep='\t', stringsAsFactors=F)
    valMat <- cbind(valMat,repVals[,3])
    valMat <- cbind(valMat,repVals[,4])
}
genenames <-repVals$gene
rld<-rlog(valMat)
print(dim(rld))

num_samples = ncol(valMat) / 2
cond_sample <- factor(rep(1:num_samples,each=2)) 
cond_ChiHum <- factor(rep(c("Human","Chimp"),times=num_samples),level=c("Human","Chimp")) #first level
cond_line<-factor(rep(str_split_fixed(sub_files, "_", 4)[,3],each=2))
sampleid<-str_split_fixed(sub_files, "_", 4)[,1]
id<-as.integer(gsub('RNA','',sampleid))
cond_batch<-c()
for (i in id){
    
    if (i<=18){
        cond_batch=c(cond_batch,1,1)
    }else{cond_batch=c(cond_batch,2,2)}
}
cond_batch<-factor(cond_batch)
celltypes<-rep(str_split_fixed(sub_files, "_", 4)[,2],each=2)
Hybs<-rep(str_split_fixed(sub_files, "_", 4)[,3],each=2)
annotation<-data.frame(allele=cond_ChiHum,hybrid=cond_line,cellType=celltypes)
data<-rld
library(pheatmap)
sabot=data[order(apply(data, 1, var), decreasing=TRUE),][1:1000,]
colnames(sabot)=paste0(rep(sub_files,each=2),rep(c('-1','-2'),num_samples))
hc<-as.dist(1 - cor(t(sabot),method='pearson'))


annoCol<-list(allele=c(Human="#cd5555", Chimp="#104e8b"),
              hybrid=c(Hyb1="#6959cd",Hyb2="#3cb371"),
              cellType=c(CM="#86B875",HP="#E495A5",PP="#C7A76C",MN="#39BEB1",SKM="#7DB0DD",RPE="#CD99D8"))
        
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots8/'
pdf(paste0(output_path,'heatmap_pearson_align2',Aligned2,'.pdf'), useDingbats=FALSE)

rownames(annotation)<-colnames(sabot)
pheatmap(as.matrix(sabot),
         show_rownames=F,
         annotation=annotation,
         show_colnames=F,
         annotation_colors = annoCol,
         clustering_distance_rows = hc,
         treeheight_row=0)
         #main='Top 1000 genes with highest variances in rlog\nGenes cluster by pearson correlation')
dev.off()

data<-rld
data = data[which(apply(data, 1, var) != 0) ,]
pca_res = prcomp(t(data), scale=TRUE, center=TRUE)
s=summary(pca_res)
df<-data.frame(pca_res$x[,1],pca_res$x[,2],celltypes,cond_ChiHum)
colnames(df)=c("PC1","PC2","cellType","Allele") 
pdf(paste0(output_path,'pca_rlog_align2',Aligned2,'.pdf'), useDingbats=FALSE)
ggplot(df, aes(x=PC1, y=PC2, shape=Allele, colour=cellType)) +
    scale_color_manual(breaks=c("CM","HP","PP","MN","SKM","RPE"),values=c("#86B875","#E495A5","#C7A76C","#39BEB1","#7DB0DD","#CD99D8")) +
    labs(x = paste0("PC1 (", round(s$importance[2]*100, 2), "%)"), y = paste0("PC2 (", round(s$importance[5]*100, 2), "%)")) +
    geom_point(size=3)+theme_bw()
dev.off()
cts<-unique(celltypes)
for (celltype in cts){
    idx<-grep(celltype,sub_files)
    idx1<-2*idx-1
    idx2<-2*idx
    data<-rld[,c(idx1,idx2)]
    celltypes1<-celltypes[c(idx1,idx2)]   
    Hybs1<-Hybs[c(idx1,idx2)]   
    cond_ChiHum1<-cond_ChiHum[c(idx1,idx2)]
    data = data[which(apply(data, 1, var) != 0) ,]
    pca_res = prcomp(t(data), scale=TRUE, center=TRUE)
    s=summary(pca_res)
    df<-data.frame(pca_res$x[,1],pca_res$x[,2],celltypes1,Hybs1,cond_ChiHum1)
    s=summary(pca_res)
    df<-data.frame(pca_res$x[,1],pca_res$x[,2],celltypes1,Hybs1,cond_ChiHum1)
    colnames(df)=c("PC1","PC2","cellType","hybridLine","Allele") 
    pdf(paste0(output_path,'pca_rlog_',celltype,'_align2',Aligned2,'.pdf'), useDingbats=FALSE)
    print(ggplot(df, aes(x=PC1, y=PC2, shape=Allele, colour=hybridLine)) +
        scale_color_manual(breaks=c("Hyb1","Hyb2"),values=c("#6959cd","#3cb371")) +
          labs(x = paste0("PC1 (", round(s$importance[2]*100, 2), "%)"), y = paste0("PC2 (", round(s$importance[5]*100, 2), "%)")) +
        geom_point(size=5)+theme_bw())
    dev.off()
}
        
    
