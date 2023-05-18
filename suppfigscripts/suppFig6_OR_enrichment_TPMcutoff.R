library("stringr")
library('ggplot2')
library(epitools)
align2 <-'human'
backup_path='/scratch/users/banwang3/celllines/data_backup'
input_path=paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',align2,'/filtered/')
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/suppplots/'
tpm<-read.csv(paste0(backup_path,input_path,'TPM.txt'),sep='\t')
cutoff=0.05 
TPMcutoffs=c(1,2,3,4,5) 
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
input_path <-paste0('/scratch/users/banwang3/celllines/analysis/DESeq2_norm/',align2,'/')
endFile<-paste0('_aligned2',align2,'_qvalue_DESeq2.txt')
result<-c()
for (ct in celltypes){
    deseq2<-read.csv(paste0(input_path,ct,'_aligned2',align2,'_qvalue_DESeq2.txt'),sep='\t')
    deseq2_chp<-read.csv(gsub("human","chimp",paste0(input_path,ct,'_aligned2',align2,'_qvalue_DESeq2.txt')),sep='\t') 
    DEgenes=merge(deseq2,deseq2_chp,by=c('gene'),suffixes=c('.hum','.chp'))
    idx=which(DEgenes$qval_mine.hum <= cutoff & DEgenes$qval_mine.chp <= cutoff & abs(DEgenes$log2FoldChange.hum-DEgenes$log2FoldChange.chp)<=1)
    DESiggenes=DEgenes[idx,]
    colidx<-c(1:length(celltypes))[-grep(ct,colnames(df))] 
    tmpdf=df[,colidx]

    print(paste0("# sig ASE genes in ",ct," : ", dim(DESiggenes)[1]))
    OR<-c()
    lower<-c()
    upper<-c()
    pvals<-c()
    for (TPMcutoff in TPMcutoffs){
        rowidx=apply(tmpdf, 1, function(x) all(x < 1)) 
        tmpdf2=df[rowidx,] 
        tmpdf3=tmpdf2[tmpdf2[,grep(ct,colnames(df))]>=TPMcutoff,]  
        print(paste0("# Tissue expressed only genes for ",ct," with TPM cutoff ",TPMcutoff, " : ", dim(tmpdf3)[1]))
        a<-dim(merge(tmpdf3,DESiggenes,by=c('gene')))[1]
        b<-dim(DESiggenes)[1]-a
        c<-dim(tmpdf3)[1]-a
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
    ORtest=data.frame(TPMcutoffs=TPMcutoffs,OR=OR,pvals=pvals,celltype=ct,lower=lower,upper=upper)
    result=rbind(result,ORtest)
}
result$plab <- ifelse(result$pvals<0.001,'***',ifelse(result$pvals<0.01,'**',ifelse(result$pvals<0.05,'*','')))
result$TPMlab<-paste0('TPM >= ',result$TPMcutoffs) 
pdf(paste0(output_path,'OR_test_tissueonlyexpressed_TPMcutoff_align2',align2,'_',SKM,'.pdf'),width=7,height=5, useDingbats=FALSE)
ggplot(data=result,aes(y=celltype)) +
    geom_point(aes(x=OR))+
    geom_segment(aes(x=lower,xend=upper,yend=celltype))+
    facet_wrap(~TPMlab,scales="free")+
    geom_vline(xintercept=1,linetype="dotted") +
    geom_text(aes(x=OR,label=plab),vjust=-0.5,fontface='bold')+
    xlab('Odds Ratio')+  
    ylab('Cell Type') +
    scale_y_discrete(limits=rev)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()

pdf(paste0(output_path,'OR_test_tissueonlyexpressed_TPMcutoff_align2',align2,'_',SKM,'_v2.pdf'),width=8,height=5, useDingbats=FALSE)
ggplot(data=result,aes(y=TPMlab))+
    geom_point(aes(x=OR))+
    geom_segment(aes(x=lower,xend=upper,yend=TPMlab))+
    facet_wrap(~celltype,scales="free")+
    geom_vline(xintercept=1,linetype="dotted") +
    geom_text(aes(x=OR,label=plab),vjust=-0.5,fontface='bold')+
    xlab('Odds Ratio')+  
    ylab('')+
    scale_y_discrete(limits=rev)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
