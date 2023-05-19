library("stringr")
#library(tidyverse)                       
library('ggplot2')
library(epitools)
align2 <-'human'
#input_path=paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',align2,'/filtered/')
#tpm<-read.csv(paste0(input_path,'TPM.txt'),sep='\t')
#print(dim(tpm))
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/'
cutoff=0.05
includeSKM=F
if (includeSKM){
    celltypes=c('CM','MN','HP','PP','RPE','SKM')
    SKM='includeSKM'
}else{
    celltypes=c('CM','MN','HP','PP','RPE')
    SKM='excludeSKM'
}


input_path <-paste0('/scratch/users/banwang3/celllines/analysis/DESeq2_norm_splitpart/',align2,'/')

endFile<-paste0('_aligned2',align2,'_qvalue_DESeq2_part_',part,'.txt')
#input_path<-paste0('/scratch/users/banwang3/celllines/analysis/DESeq2_norm/',align2,'/')
#endFile<-paste0('_aligned2',align2,'_qvalue_DESeq2.txt')    
#input_files<-list.files(input_path)
#input_files<-input_files[grep("_vs_",input_files)]
parts=c(1,2)
result<-c()
for (part in parts){
    for (part1 in parts){
OR<-c()
lower<-c()
upper<-c()
pvals<-c()

for (ct in celltypes){
    deseq2<-read.csv(paste0(input_path,ct,'_aligned2',align2,'_qvalue_DESeq2_part_',part,'.txt'),sep='\t')
    deseq2_chp<-read.csv(gsub("human","chimp",paste0(input_path,ct,'_aligned2',align2,'_qvalue_DESeq2_part_',part,'.txt')),sep='\t') 
    DEgenes=merge(deseq2,deseq2_chp,by=c('gene'),suffixes=c('.hum','.chp'))
    idx=which(DEgenes$qval_mine.hum <= cutoff & DEgenes$qval_mine.chp <= cutoff & abs(DEgenes$log2FoldChange.hum-DEgenes$log2FoldChange.chp)<=1)
    DESiggenes=DEgenes[idx,]
    print(paste0("# total genes for ",ct," :", dim(DEgenes)[1]))
    print(paste0("# DE genes for ",ct," :", dim(DESiggenes)[1]))  
    TSmerge=c()
    #part1 = parts[parts!=part] 
    endFile<-paste0('_aligned2',align2,'_qvalue_DESeq2_part_',part1,'.txt')
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
        a<-dim(merge(TSmerge,DESiggenes,by=c('gene')))[1]
        b<-dim(DESiggenes)[1]-a
        c<-dim(TSmerge)[1]-a
        d<-dim(DEgenes)[1]-a-b-c
        dat=data.frame('TS'=c(a,c),'nTS'=c(b,d),row.names = c("ASE", "nASE"),stringsAsFactors = FALSE)
        ORtable<-matrix(c(a,c,b,d),nrow = 2, ncol = 2) 
        result1=oddsratio.wald(ORtable)
        OR<-c(OR,result1$measure[2,1])
        lower<-c(lower,result1$measure[2,2])
        upper<-c(upper,result1$measure[2,3])
        pvals<-c(pvals,result1$p.value[2,3])
        print(dat)
    }
    ORtest=data.frame(celltype=celltypes,OR=OR,pvals=pvals,lower=lower,upper=upper,DEpart=part,TSpart=part1)
    result=rbind(result,ORtest)
}
}
result$plab <- ifelse(result$pvals<0.001,'***',ifelse(result$pvals<0.01,'**',ifelse(result$pvals<0.05,'*','')))
result$partlabel<-factor(paste0('DEpart-',result$DEpart,' & ','TSpart-',result$TSpart))
#result$pdisplay<-ifelse(result$pvals<cutoff,signif(result$pvals,digits=1),'not sig')    
#result1<-result[result$TPMcutoffs==5,]
pdf(paste0(output_path,'OR_test_splitparts_',SKM,'.pdf'), useDingbats=FALSE)
ggplot(data=result,aes(y=celltype)) +
    geom_point(aes(x=OR))+
    geom_segment(aes(x=lower,xend=upper,yend=celltype))+
    ylab('Cell Type') +
    xlab('Odds Ratio')+
    facet_wrap(~partlabel) +
    geom_vline(xintercept=1,linetype="dotted") +
    geom_text(aes(x=OR,label=plab),vjust=-0.5,fontface='bold')+
    scale_y_discrete(limits=rev)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
