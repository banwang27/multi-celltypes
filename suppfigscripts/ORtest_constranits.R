library("stringr")
library('ggplot2')
library(epitools)
align2 <-'human'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/suppplots/'
includeSKM=T
if (includeSKM){
    celltypes=c('CM','MN','HP','PP','RPE','SKM')
    SKM='includeSKM'
}else{
    celltypes=c('CM','MN','HP','PP','RPE')
    SKM='excludeSKM'
}
df = read.csv('constrains/Haploinsufficiency_Gene_Predictions.csv',sep=',')
df=df[,c('Gene','pHaplo','gnomAD.pLI')]
colnames(df)=c('gene','pHI','pLI')
row.has.na <- apply(df, 1, function(x){any(is.na(x))})
df <- df[!row.has.na,]
allfiles = list.files('constrains/')
input_path <-paste0('/scratch/users/banwang3/celllines/analysis/DESeq2_norm/',align2,'/')
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
    ASEvar = read.csv(paste0('constrains/',allfiles[grep(ct,allfiles)]),sep='\t')
    ASEvar = ASEvar[,c('gene','Variance.of.ASE.Dist')]
    colnames(ASEvar)=c('gene','ASEvar')
    df1 = merge(df,ASEvar,all.x=T,all.y=T)
    for (keyword in c('pHI','pLI','ASEvar')){
        row.has.na<-which(is.na(df1[,keyword])) 
        df2<-df1[-row.has.na,]
        df2 = df2[order(df2[,keyword]),] 
        df2 = df2[df2$gene %in% DEgenes$gene,] 
        steps=floor(nrow(df2)/5)
        OR<-c()
        lower<-c()
        upper<-c()
        pvals<-c()
        for (bin in seq(1,5)){
            if (bin==5){
                df3 = df2[(4*steps+1):nrow(df2),]
            }else{
                df3 = df2[((bin-1)*steps+1):(bin*steps),]
            }
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
        }
        ORtest=data.frame(celltype=ct,OR=OR,pvals=pvals,lower=lower,upper=upper,expressionbin=seq(1,5),keyword=keyword)
        result=rbind(result,ORtest)
    }
}
result$plab <- ifelse(result$pvals<0.001,'***',ifelse(result$pvals<0.01,'**',ifelse(result$pvals<0.05,'*','')))
result$binlabshort<-NA
result[result$expressionbin==1,'binlabshort']<-'0-20%'
result[result$expressionbin==2,'binlabshort']<-'20-40%'
result[result$expressionbin==3,'binlabshort']<-'40-60%'
result[result$expressionbin==4,'binlabshort']<-'60-80%'
result[result$expressionbin==5,'binlabshort']<-'80-100%'
result$binlabshort=factor(result$binlabshort,levels=c('0-20%','20-40%','40-60%','60-80%','80-100%'))

for (keyword in c('pHI','pLI','ASEvar')){
    pdf(paste0(output_path,'OR_test_',keyword,'_',SKM,'_.pdf'), useDingbats=FALSE,width=8,height=5)
    result1 = result[result$keyword==keyword,]
    p<-ggplot(data=result1,aes(y=binlabshort))+
    geom_point(aes(x=OR))+geom_segment(aes(x=lower,xend=upper,yend=binlabshort))+
    ylab(paste0('Genes with top% value in ',keyword))+
    xlab('Odds Ratio')+
    facet_wrap(~celltype,scales="free")+
    geom_vline(xintercept=1,linetype="dotted") +
    geom_text(aes(x=OR,label=plab),vjust=-0.5,fontface='bold')+
    scale_y_discrete(limits=rev)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    print(p)
    dev.off()
}
