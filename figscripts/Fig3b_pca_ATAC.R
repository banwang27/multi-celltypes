library(DESeq2)
library(stringr)
library(ggplot2)
library(pheatmap)
align2='human'
reffed='humreffed'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/'
dir_folder=paste0('/scratch/users/banwang3/celllines/ATAC/All_Peaks_counts/')
filename='_ATAC_Filtered_all_peaks.txt'
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
celltypes=c('CM','HP','MN','PP','SKM')
celltype=celltypes[1]
filepath=paste0(dir_folder,celltype,filename)
df1=read.csv(filepath,sep='\t',header=T)  
idx = grep(reffed,colnames(df1))   
df = df1[,c(1,idx)]
for (celltype in celltypes[2:length(celltypes)]){
    filepath=paste0(dir_folder,celltype,filename)
    df1=read.csv(filepath,sep='\t',header=T)
    idx = grep(reffed,colnames(df1))
    df = merge(df,df1[,c(1,idx)])
}
ASE_table=as.matrix(df[,2:ncol(df)])
valMat = ASE_table
rld<-rlog(valMat)
print(dim(rld))
num_samples = ncol(valMat) / 2
sub_files = colnames(ASE_table)
celltypes<-str_split_fixed(sub_files, "_", 7)[,2]
Hybs<-str_split_fixed(sub_files, "_", 7)[,3]
cond_ChiHum<-str_split_fixed(sub_files, "_", 7)[,4] 
cond_ChiHum<-factor(cond_ChiHum,level=c("Human","Chimp"))
annotation<-data.frame(allele=cond_ChiHum,hybrid=Hybs,cellType=celltypes)
data<-rld
sabot=data[order(apply(data, 1, var), decreasing=TRUE),][1:1000,]
colnames(sabot)=sub_files
hc<-as.dist(1 - cor(t(sabot),method='pearson'))


annoCol<-list(allele=c(Human="indianred1", Chimp="dodgerblue1"),
              hybrid=c(Hyb1="#6959cd",Hyb2="#3cb371"),
              cellType=c(CM="#86B875",HP="#E495A5",PP="#C7A76C",MN="#39BEB1",SKM="#7DB0DD"))
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/'
pdf(paste0(output_path,'heatmap_pearson_rlog_ATAC_align2',align2,'.pdf'), useDingbats=FALSE)

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

data<-rld[,c(1:ncol(valMat))]
data<-rld
data = data[which(apply(data, 1, var) != 0) ,]
pca_res = prcomp(t(data), scale=TRUE, center=TRUE)
s=summary(pca_res)
df<-data.frame(pca_res$x[,1],pca_res$x[,2],celltypes,cond_ChiHum)
colnames(df)=c("PC1","PC2","cellType","Allele") 
pdf(paste0(output_path,'pca_rlog_ATAC_align2',align2,'.pdf'), width=4.5,height=3.5,useDingbats=FALSE)
ggplot(df, aes(x=PC1, y=PC2, shape=Allele, colour=cellType)) +
    scale_color_manual(breaks=c("CM","HP","PP","MN","SKM"),values=c("#86B875","#E495A5","#C7A76C","#39BEB1","#7DB0DD")) +
    labs(x = paste0("PCA 1 (", round(s$importance[2]*100, 2), "%)"), y = paste0("PCA 2 (", round(s$importance[5]*100, 2), "%)")) +
    geom_point(size=3)+theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 
dev.off()
cts<-unique(celltypes)
for (celltype in cts){
    idx<-grep(celltype,sub_files)
    data<-rld[,idx]
    celltypes1<-celltypes[idx]   
    Hybs1<-Hybs[idx]   
    cond_ChiHum1<-cond_ChiHum[idx]
    data = data[which(apply(data, 1, var) != 0) ,]
    pca_res = prcomp(t(data), scale=TRUE, center=TRUE)
    s=summary(pca_res)
    df<-data.frame(pca_res$x[,1],pca_res$x[,2],celltypes1,Hybs1,cond_ChiHum1)
    s=summary(pca_res)
    df<-data.frame(pca_res$x[,1],pca_res$x[,2],celltypes1,Hybs1,cond_ChiHum1)
    colnames(df)=c("PC1","PC2","cellType","hybridLine","Allele") 
    pdf(paste0(output_path,'pca_rlog_',celltype,'_ATAC_align2',align2,'.pdf'), width=4.5, height=3.5,useDingbats=FALSE)
    print(ggplot(df, aes(x=PC1, y=PC2, shape=Allele, colour=hybridLine)) +
        scale_color_manual(breaks=c("Hyb1","Hyb2"),values=c("#6959cd","#3cb371")) +
          labs(x = paste0("PCA 1 (", round(s$importance[2]*100, 2), "%)"), y = paste0("PCA 2 (", round(s$importance[5]*100, 2), "%)")) +
            geom_point(size=3)+theme_bw()+
            theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
    dev.off()
}
