library("stringr")
library('ggplot2')
library(epitools)
align2 <-'human'
backup_path='/scratch/users/banwang3/celllines/data_backup'
input_path=paste0(backup_path,'/scratch/users/banwang3/celllines/analysis/ASE_norm/',align2,'/filtered/')
tpm<-read.csv(paste0(input_path,'TPM.txt'),sep='\t')
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/'
includeSKM=T
if (includeSKM){
    celltypes=c('CM','MN','HP','PP','RPE','SKM')
    SKM='includeSKM'
    totaltpm<-tpm[,grep('total',colnames(tpm))] 
}else{
    celltypes=c('CM','MN','HP','PP','RPE')
    SKM='excludeSKM'
    totaltpm<-tpm[,grep('total',colnames(tpm))] 
    totaltpm<-totaltpm[,-grep('SKM',colnames(totaltpm))]
}
df=c()
for (ct in celltypes){
    totaltpm1<-totaltpm[,grep(ct,colnames(totaltpm))] 
    df = cbind(df,rowMeans(totaltpm1)) 
}
df=as.data.frame(df)
colnames(df)=celltypes
df$gene=tpm$gene

input_path <-paste0(backup_path,'/scratch/users/banwang3/celllines/analysis/DESeq2_norm/',align2,'/')
endFile<-paste0('_aligned2',align2,'_qvalue_DESeq2.txt')
cutoff=0.05
result<-c()
for (ct in celltypes){
        deseq2<-read.csv(paste0(input_path,ct,'_aligned2',align2,'_qvalue_DESeq2.txt'),sep='\t')
        deseq2_chp<-read.csv(gsub("human","chimp",paste0(input_path,ct,'_aligned2',align2,'_qvalue_DESeq2.txt')),sep='\t') 
        DEgenes=merge(deseq2,deseq2_chp,by=c('gene'),suffixes=c('.hum','.chp'))
        idx=which(DEgenes$qval_mine.hum <= cutoff & DEgenes$qval_mine.chp <= cutoff & abs(DEgenes$log2FoldChange.hum-DEgenes$log2FoldChange.chp)<=1)
        DESiggenes=DEgenes[idx,]
        print(paste0("# total genes for ",ct," :", dim(DEgenes)[1]))
        print(paste0("# DE genes for ",ct," :", dim(DESiggenes)[1]))  
        TSmerge=c()
        endFile<-paste0('_aligned2',align2,'_qvalue_DESeq2.txt')
        for (ct1 in celltypes[celltypes!=ct]){
            TSdeseq2<-read.csv(paste0(input_path,ct,'_vs_',ct1,endFile),sep='\t')
            TSdeseq2_chp<-read.csv(gsub("human","chimp",paste0(input_path,ct,'_vs_',ct1,endFile)),sep='\t') 
            TSgenes=merge(TSdeseq2,TSdeseq2_chp,by=c('gene'),suffixes=c('.hum','.chp'))
            idx=which(TSgenes$qval_mine.hum <= cutoff & TSgenes$qval_mine.chp <= cutoff & abs(TSgenes$log2FoldChange.hum-TSgenes$log2FoldChange.chp)<=1)  
            TSgenes=TSgenes[idx,] 
            TSgenes[,ct1]=1
            if (length(TSmerge)==0){
                TSmerge=TSgenes[,c('gene',ct1)]
            }else{
                TSmerge=merge(TSmerge,TSgenes[,c('gene',ct1)],by=c('gene'))
            }
        }   
        print(paste0("# TS genes for ",ct," :", dim(TSmerge)[1]))   
        df1 = df[order(df[,ct]),] 
        df3 = df1[df1$gene %in% DEgenes$gene,] 
        OR<-c()
        lower<-c()
        upper<-c()
        pvals<-c()
        a<-dim(merge(TSmerge[TSmerge$gene%in%df3$gene,],DESiggenes[DESiggenes$gene%in%df3$gene,],by=c('gene')))[1]
        b<-dim(DESiggenes[DESiggenes$gene%in%df3$gene,])[1]-a
        c<-dim(TSmerge[TSmerge$gene%in%df3$gene,])[1]-a
        d<-dim(DEgenes[DEgenes$gene%in%df3$gene,])[1]-a-b-c
        dat=data.frame('TS'=c(a,c),'nTS'=c(b,d),row.names = c("ASE", "nASE"),stringsAsFactors = FALSE)
        ORtable<-matrix(c(a,c,b,d),nrow = 2, ncol = 2) 
        result1=oddsratio.wald(ORtable)
        OR<-c(OR,result1$measure[2,1])
        lower<-c(lower,result1$measure[2,2])
        upper<-c(upper,result1$measure[2,3])
        pvals<-c(pvals,result1$p.value[2,3])
        print(dat)
    ORtest=data.frame(celltype=ct,OR=OR,pvals=pvals,lower=lower,upper=upper)
    result=rbind(result,ORtest)
}
result$plab <- ifelse(result$pvals<0.001,'***',ifelse(result$pvals<0.01,'**',ifelse(result$pvals<0.05,'*','')))
result$celltypefull=ifelse(result$celltype=='CM','Cardiomyocyte',ifelse(result$celltype=='HP','Hepatocyte',ifelse(result$celltype=='MN','Motor\nNeuron',ifelse(result$celltype=='PP','Pancreatic\nProgenitor',ifelse(result$celltype=='RPE','Retinal\nPigment\nEpithelium','Skeletal\nMyocyte')))))
pdf(paste0(output_path,'OR_test_',SKM,'.pdf'), useDingbats=FALSE,width=4.5,height=4)
ggplot(data=result,aes(y=celltypefull)) +
    geom_point(aes(x=OR))+
    geom_segment(aes(x=lower,xend=upper,yend=celltypefull))+
    ylab('') +
    xlab('Odds Ratio')+
    geom_vline(xintercept=1,linetype="dotted") +
    geom_text(aes(x=OR,label=plab),vjust=-0.5,fontface='bold')+
    scale_y_discrete(limits=rev)+
    xlim(0.5,2.5)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()

