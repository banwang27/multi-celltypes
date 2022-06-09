library(stringr)
library(ggplot2)
library(corrplot)
library(pheatmap)
align2='human'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots3/'
dir_folder=paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',align2,'/')
samples=c('RNA3_CM_Hyb1','RNA4_CM_Hyb2','RNA7_SKM_Hyb1','RNA8_SKM_Hyb2','RNA11_HP_Hyb1','RNA12_HP_Hyb2','RNA15_PP_Hyb1','RNA16_PP_Hyb2','RNA17_CM_Hyb1','RNA18_MN_Hyb2')
samples2 = c('RNA19_CM_Hyb2','RNA20_SKM_Hyb1','RNA21_SKM_Hyb2','RNA22_HP_Hyb1','RNA23_HP_Hyb2','RNA24_PP_Hyb1','RNA25_PP_Hyb2','RNA26_RPE_Hyb1','RNA27_RPE_Hyb1','RNA28_RPE_Hyb2','RNA29_RPE_Hyb2','RNA30_MN_Hyb1','RNA31_MN_Hyb2')
samples_all = c(samples,samples2)
filename='_ase_normalized.txt'
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
##Plot PCA on human and chimp allele
ASE_table=c()
total=c()
for (sample1 in samples_all){
    filepath=paste0(dir_folder,sample1,filename)
    df=read.csv(filepath,sep='\t')
    df_total=sum(colSums(df[,3:6]))
    ASE_table=cbind(ASE_table,df$ref_counts,df$alt_counts)
    total=c(total,df_total,df_total) 
}
#Normalize counts to CPM
ASE_table= sweep(as.matrix(ASE_table), 2, as.vector(total/1000000), `/`)
#Filter out low expressed genes
ASE_table_1 = ASE_table[which(rowSums(ASE_table)>1),]
#Perform PCA
pca_res = prcomp(t(ASE_table_1), scale=TRUE, center=TRUE)
s=summary(pca_res)
#Define metadata
#sampleID
tmp <- rep(samples_all,each=2)
#cellType
my_group<-as.factor(str_split(tmp, '_',simplify = TRUE)[,2])
#Species
species<-rep(c('hum','chp'),length(tmp)/2)
#Plot PCA
df<-data.frame(pca_res$x[,1],pca_res$x[,2],tmp,my_group,species)
colnames(df)=c("PC1","PC2","sampleID","cellType","Allele") 
pdf(paste0(output_path,'pca_both_allele_align2',align2,'.pdf'), useDingbats=FALSE)
ggplot(df, aes(x=PC1, y=PC2, shape=Allele, colour=cellType)) +
  labs(x = paste0("PCA 1 (", round(s$importance[2]*100, 2), "%)"), y = paste0("PCA 1 (", round(s$importance[5]*100, 2), "%)"), title = "PCA on both allele counts") +
  geom_point(size=3)+theme_bw()
dev.off()


##Plot PCA on human or chimp allele in one cell type
celltypes=unique(str_split(samples_all, '_',simplify = TRUE)[,2])
for (celltype in celltypes){
  ASE_table=c()
  celltypesamples=samples_all[grepl(celltype, samples_all)]
  total=c()
  for (sample1 in celltypesamples){
    filepath=paste0(dir_folder,sample1,filename)
    df=read.csv(filepath,sep='\t')
    df_total=sum(colSums(df[,3:6]))
    ASE_table=cbind(ASE_table,df$ref_counts,df$alt_counts)  
    total=c(total,df_total,df_total) 
  }
  ASE_table= sweep(as.matrix(ASE_table), 2, as.vector(total/1000000), `/`)
  ASE_table_1 = ASE_table[which(rowSums(ASE_table)>1),]
  pca_res = prcomp(t(ASE_table_1), scale=TRUE, center=TRUE)
  s=summary(pca_res)  
  tmp <- rep(celltypesamples,each=2)
  my_group<-str_split(tmp, '_',simplify = TRUE)[,2]
  Hyb<-as.factor(str_split(tmp, '_',simplify = TRUE)[,3])
  species<-rep(c('hum','chp'),length(tmp)/2)
  df<-data.frame(pca_res$x[,1],pca_res$x[,2],tmp,Hyb,species)
  colnames(df)=c("PC1","PC2","sampleID","hybridLine","Allele") 
  pdf(paste0(output_path,'pca_both_allele_',celltype,'_align2',align2,'.pdf'), useDingbats=FALSE)
  print(ggplot(df, aes(x=PC1, y=PC2, shape=Allele, colour=hybridLine)) +
          labs(x = paste0("PCA 1 (", round(s$importance[2]*100, 2), "%)"), y = paste0("PCA 1 (", round(s$importance[5]*100, 2), "%)"), title = paste0("PCA on both allele counts ",celltype)) +
          geom_point(size=3)+theme_bw())
  dev.off()
}

##Plot PCA on ASE
ASE_table=c()
ASE_table1=c()
for (sample1 in samples_all){
  filepath=paste0(dir_folder,sample1,filename)
  df=read.csv(filepath,sep='\t')
  df$ASE=NA
  idx=which(df$ref_counts>=1 & df$alt_counts>=1)
  df[idx,'ASE']=log2(df[idx,'ref_counts']/df[idx,'alt_counts'])
  ASE_table=cbind(ASE_table,df$ASE)
  ASE_table1=cbind(ASE_table1,df$ref_counts,df$alt_counts)
}
ASE_table1= sweep(as.matrix(ASE_table1), 2, as.vector(colSums(ASE_table1)/1000000), `/`)
idx=which(rowSums(ASE_table1)<=1)
ASE_table[idx,]=NA
row.has.na <- apply(ASE_table, 1, function(x){any(is.na(x))})
data <- ASE_table[!row.has.na,]
tmp<-samples_all
Hyb<-as.factor(str_split(tmp, '_',simplify = TRUE)[,3])
my_group<-as.factor(str_split(tmp, '_',simplify = TRUE)[,2])
#Batch
Batch<-c(rep("prep1",10),rep("prep2",13))
a = data.frame(t(data))
pca_res<-prcomp(a,scale=TRUE) 
s<-summary(pca_res)
df<-data.frame(pca_res$x[,1],pca_res$x[,2],tmp,my_group,Hyb)
colnames(df)=c("PC1","PC2","sampleID","cellType","hybridLine") 
pdf(paste0(output_path,'pca_ASE_align2',align2,'.pdf'), useDingbats=FALSE)
ggplot(df, aes(x=PC1, y=PC2, shape=hybridLine, colour=cellType)) +
  labs(x = paste0("PCA 1 (", round(s$importance[2]*100, 2), "%)"), y = paste0("PCA 1 (", round(s$importance[5]*100, 2), "%)"), title = "PCA on ASE log2(hum/chp)") +
  geom_point(size=3)+theme_bw()
dev.off()

#Plot heatmap
abot<-as.matrix(data)
annotation<-data.frame(cellType=my_group,hybrid=Hyb,batch=Batch)
print(dim(abot))
sabot = abot[order(apply(abot, 1, var), decreasing=TRUE),][1:1000,]
colnames(sabot)=samples_all
hc<-as.dist(1 - cor(t(sabot),method='pearson'))
pdf(paste0(output_path,'heatmap_pearson_align2',align2,'.pdf'), useDingbats=FALSE)
rownames(annotation)<-samples_all
pheatmap(as.matrix(sabot),show_rownames=F,clustering_distance_rows = hc,main='Top 1000 genes with highest variances in log2(hum/chp) values\nGenes cluster by pearson correlation',annotation=annotation)
dev.off()

hc<-as.dist(1 - cor(t(sabot),method='spearman'))
pdf(paste0(output_path,'heatmap_spearman_align2',align2,'.pdf'), useDingbats=FALSE)
rownames(annotation)<-colnames(sabot)
pheatmap(as.matrix(sabot),show_rownames=F,clustering_distance_rows = hc,main='Top 1000 genes with highest variances in log2(hum/chp) values\nGenes cluster by spearman correlation',annotation=annotation)
dev.off()

#Plot correlation plot
data <- as.matrix(data)
M = cor(data,method='spearman')
row.names(M)=samples_all
colnames(M)=samples_all 
pdf(paste0(output_path,'corplot_ASE_spearman_align2',align2,'.pdf'), useDingbats=FALSE)
corrplot(M, method = 'color',type='upper',tl.cex=0.4,addCoef.col = 'black',number.cex = 0.4,is.corr = FALSE,order = "hclust",tl.col='black',col= colorRampPalette(c("white","pink","red"))(10),col.lim = c(floor_dec(min(M),1),1))
dev.off()

M = cor(data,method='pearson')
row.names(M)=samples_all
colnames(M)=samples_all 
pdf(paste0(output_path,'corplot_ASE_pearson_align2',align2,'.pdf'), useDingbats=FALSE)
corrplot(M, method = 'color',type='upper',tl.cex=0.4,addCoef.col = 'black',number.cex = 0.4,is.corr = FALSE,order = "hclust",tl.col='black',col= colorRampPalette(c("white","pink","red"))(10),col.lim = c(floor_dec(min(M),1),1))
dev.off()

##Plot PCA on all counts (ase+noase)
ASE_table=c()
for (sample1 in samples_all){
  filepath=paste0(dir_folder,sample1,filename)
  df=read.csv(filepath,sep='\t')
  df$ASE=df$ref_counts+df$alt_counts+df$no_ase_counts+df$ambig_ase_counts
  ASE_table=cbind(ASE_table,df$ASE)
}
ASE_table= sweep(as.matrix(ASE_table), 2, as.vector(colSums(ASE_table)/1000000), `/`)
ASE_table_1 = ASE_table[which(rowSums(ASE_table)>1),]
pca_res = prcomp(t(ASE_table_1), scale=TRUE, center=TRUE)
s=summary(pca_res)
tmp <- samples_all
Hyb<-as.factor(str_split(tmp, '_',simplify = TRUE)[,3])
my_group<-as.factor(str_split(tmp, '_',simplify = TRUE)[,2])
#Batch
Batch<-c(rep("prep1",10),rep("prep2",13))
df<-data.frame(pca_res$x[,1],pca_res$x[,2],tmp,my_group,Hyb)
colnames(df)=c("PC1","PC2","sampleID","cellType","hybridLine") 
pdf(paste0(output_path,'pca_all_counts_align2',align2,'.pdf'), useDingbats=FALSE)
ggplot(df, aes(x=PC1, y=PC2, shape=hybridLine, colour=cellType)) +
  labs(x = paste0("PCA 1 (", round(s$importance[2]*100, 2), "%)"), y = paste0("PCA 1 (", round(s$importance[5]*100, 2), "%)"), title = "PCA on all counts (ASE+ no ASE)") +
  geom_point(size=3)+ theme_bw()
dev.off()

#Plot heatmap
abot<-as.matrix(ASE_table_1)
annotation<-data.frame(cellType=my_group,hybrid=Hyb,batch=Batch)
print(dim(abot))
sabot = abot[order(apply(abot, 1, var), decreasing=TRUE),][1:1000,]
colnames(sabot)=samples_all
hc<-as.dist(1 - cor(t(sabot),method='pearson'))
pdf(paste0(output_path,'heatmap_pearson_align2',align2,'.pdf'), useDingbats=FALSE)
rownames(annotation)<-samples_all
pheatmap(as.matrix(sabot),show_rownames=F,clustering_distance_rows = hc,main='Top 1000 genes with highest variances in log2(hum/chp) values\nGenes cluster by pearson correlation',annotation=annotation)
dev.off()

hc<-as.dist(1 - cor(t(sabot),method='spearman'))
pdf(paste0(output_path,'heatmap_spearman_align2',align2,'.pdf'), useDingbats=FALSE)
rownames(annotation)<-colnames(sabot)
pheatmap(as.matrix(sabot),show_rownames=F,clustering_distance_rows = hc,main='Top 1000 genes with highest variances in log2(hum/chp) values\nGenes cluster by spearman correlation',annotation=annotation)
dev.off()

#Plot correlation plot
data <- as.matrix(ASE_table_1)
M = cor(data,method='spearman')
row.names(M)=samples_all
colnames(M)=samples_all 
pdf(paste0(output_path,'corplot_allcounts_spearman_align2',align2,'.pdf'), useDingbats=FALSE)
corrplot(M, method = 'color',type='upper',tl.cex=0.4,addCoef.col = 'black',number.cex = 0.4,is.corr = FALSE,order = "hclust",tl.col='black',col= colorRampPalette(c("white","pink","red"))(10),col.lim = c(floor_dec(min(M),1),1))
dev.off()

M = cor(data,method='pearson')
row.names(M)=samples_all
colnames(M)=samples_all 
pdf(paste0(output_path,'corplot_allcounts_pearson_align2',align2,'.pdf'), useDingbats=FALSE)
corrplot(M, method = 'color',type='upper',tl.cex=0.4,addCoef.col = 'black',number.cex = 0.4,is.corr = FALSE,order = "hclust",tl.col='black',col= colorRampPalette(c("white","pink","red"))(10),col.lim = c(floor_dec(min(M),1),1))
dev.off()