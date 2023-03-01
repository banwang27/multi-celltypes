library(stringr)
library(ggplot2)
library(dplyr)
library(epitools)
celltype='MN'
celltypes=c('CM','MN','HP','PP','SKM')
reffed='humreffed'
align2='human'
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/suppplots/'
lfccutoff=0.5
final=c()
df=read.csv('/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/alexdata/ATAC_All_Peaks_Gnomad_Constraint_Summary.txt',sep='\t')
df = df[,1:2]
colnames(df)=c('peak_id','maxZ')
row.has.na <- apply(df, 1, function(x){any(is.na(x))})
constrains <- df[!row.has.na,]
buildORtable<-function(df,keyword,celltype,cutofftype,cutoff){
    df = df[order(df[,keyword]),]
    df = df[!duplicated(df$peak_id),]
    if (cutofftype=='lfc'){
        idx = which(abs(df$lfc)>=cutoff)
    }else if (cutofftype=='pval'){
        idx = which(df$pval<=cutoff)
    }else{
        idx = which(df$FDR<=cutoff)
    }
    dCA = df[idx,]
    celltypes1=celltypes[celltypes!=celltype] 
    idx_onlyct=grepl(celltype,df$peak_id)
    idx_nootherct=!grepl(paste(celltypes1, collapse="|"), df$peak_id)
    df_ctspecific = df[idx_onlyct & idx_nootherct,]
    OR<-c()
    lower<-c()
    upper<-c()
    pvals<-c()
    bins<-c()
    steps=floor(nrow(df)/5)
    for (bin in seq(1,5)){
         if (bin==5){
             df3 = df[(4*steps+1):nrow(df),]
         }else{
             df3 = df[((bin-1)*steps+1):(bin*steps),]
         }
         a<-dim(merge(df_ctspecific[df_ctspecific$peak_id %in% df3$peak_id,],dCA[dCA$peak_id %in% df3$peak_id,],by=c('peak_id')))[1]
         b<-dim(dCA[dCA$peak_id %in% df3$peak_id,])[1]-a
         c<-dim(df_ctspecific[df_ctspecific$peak_id %in% df3$peak_id,])[1]-a
         d<-dim(df[df$peak_id %in% df3$peak_id,])[1]-a-b-c
         ORtable<-matrix(c(a,c,b,d),nrow = 2, ncol = 2) 
         if (d>0){
            result1=oddsratio.wald(ORtable)
            OR<-c(OR,result1$measure[2,1])
            lower<-c(lower,result1$measure[2,2])
            upper<-c(upper,result1$measure[2,3])
            pvals<-c(pvals,result1$p.value[2,3])
            bins<-c(bins,bin)
        }
    }
    if (length(OR)>0){
        ORtest=data.frame(celltype=celltype,OR=OR,pvals=pvals,lower=lower,upper=upper,bin=bins,keyword=keyword,cutofftype=cutofftype,cutoff=cutoff)
    }else{
        ORtest=NA
    }
    return(ORtest)
}
for (cutofftype in c('lfc','pval','FDR')){
    if (cutofftype=='pval'){
        cutoff=0.05
    }else if(cutofftype=='lfc'){
        cutoff=0.5
    }else{
        cutoff=0.05
    }
    result=c()
    for (celltype in celltypes){
        dfopen=read.csv(paste0(input_path,'Filtered/',celltype,'_ATAC_Filtered_All_Peaks_fixed.txt'),sep='\t')
        dfopen = dfopen[,c('Peak.Humreffed','l2fc_Humr','Human.binom_pval','Human.FDR')]
        colnames(dfopen)=c('peak_id','lfc','pval','FDR')
        dfopen$celltype=celltype
        df1 = merge(dfopen,constrains)
        keyword='maxZ'
        df2 = df1[order(df1[,keyword]),] 
    
    result=rbind(result,buildORtable(df2,keyword,celltype,cutofftype,cutoff))
    }
    result$binlabshort<-NA
    result[result$bin==1,'binlabshort']<-'0-20%'
    result[result$bin==2,'binlabshort']<-'20-40%'
    result[result$bin==3,'binlabshort']<-'40-60%'
    result[result$bin==4,'binlabshort']<-'60-80%'
    result[result$bin==5,'binlabshort']<-'80-100%'
    result$binlabshort=factor(result$binlabshort,levels=c('0-20%','20-40%','40-60%','60-80%','80-100%'))
    result$plab <- ifelse(result$pvals<0.001,'***',ifelse(result$pvals<0.01,'**',ifelse(result$pvals<0.05,'*','')))
    pdf(paste0(output_path,'ATAC_OR_test_',keyword,'_cutoff_',cutofftype,'_',cutoff,'.pdf'), useDingbats=FALSE,width=8,height=5)
    p<-ggplot(data=result,aes(y=binlabshort))+    
        geom_point(aes(x=OR))+geom_segment(aes(x=lower,xend=upper,yend=binlabshort))+
        ylab(paste0('Peaks with top% value in ',keyword))+
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
