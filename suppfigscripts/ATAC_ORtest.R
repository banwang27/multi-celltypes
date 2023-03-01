library(stringr)
library(ggplot2)
library(dplyr)
library(epitools)
celltypes=c('CM','MN','HP','PP','SKM')
reffed='humreffed'
align2='human'
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/suppplots/'
result=c()
pvalcutoffs=c(0.1,0.05,0.01)
FDRcutoffs=c(0.1,0.05,0.01)
lfccutoffs=c(0.5,1,1.5)
for (celltype in celltypes){
    dfopen = read.csv(paste0(input_path,'Filtered/',celltype,'_ATAC_Filtered_All_Peaks_fixed.txt'),sep='\t')
    dfopen = dfopen[,c('Peak.Humreffed','l2fc_Humr','Human.binom_pval','Human.FDR')]    
    colnames(dfopen)=c('peak_id','lfc','pval','FDR')
    dfopen$celltype=celltype
    for (i in seq(1,3)){
        idx1 = which(dfopen$pval<=pvalcutoffs[i])
        idx2= which(dfopen$FDR<=FDRcutoffs[i])
        idx3 = which(abs(dfopen$lfc)>=lfccutoffs[i])
        dCA_pval=dfopen[idx1,]
        dCA_FDR=dfopen[idx2,]
        dCA_lfc=dfopen[idx3,]
        celltypes1=celltypes[celltypes!=celltype]
        idx_onlyct=grepl(celltype,dfopen$peak_id)
        idx_nootherct=!grepl(paste(celltypes1, collapse="|"), dfopen$peak_id)
        dfopen_ctspecific = dfopen[idx_onlyct & idx_nootherct,]
        OR<-c()
        lower<-c()
        upper<-c()
        pvals<-c()
        keywords<-c()
        cutoffs<-c()
        for (keyword in c('pval','FDR','lfc')){
            if (keyword=='pval'){
                dCA=dCA_pval
                cutoff=pvalcutoffs[i]
            }else if (keyword=='FDR'){
                dCA=dCA_FDR
                cutoff=FDRcutoffs[i]
            }else{
                dCA=dCA_lfc
                cutoff=lfccutoffs[i]
            }
            a<-dim(merge(dfopen_ctspecific,dCA,by=c('peak_id')))[1]
            b<-dim(dCA)[1]-a
            c<-dim(dfopen_ctspecific)[1]-a
            d<-dim(dfopen)[1]-a-b-c
            dat=data.frame('TS'=c(a,c),'nTS'=c(b,d),row.names = c("ASCA", "nASCA"),stringsAsFactors = FALSE)
            ORtable<-matrix(c(a,c,b,d),nrow = 2, ncol = 2) 
            result1=oddsratio.wald(ORtable)
            OR<-c(OR,result1$measure[2,1])
            lower<-c(lower,result1$measure[2,2])
            upper<-c(upper,result1$measure[2,3])
            pvals<-c(pvals,result1$p.value[2,3])
            keywords<-c(keywords,keyword)
            cutoffs<-c(cutoffs,cutoff)
            print(keyword)
            print(cutoff)
            print(dat)
        }
        ORtest=data.frame(celltype=celltype,OR=OR,pvals=pvals,lower=lower,upper=upper,cutoff=cutoffs,keyword=keywords)
        result=rbind(result,ORtest)
    }
}
result$plab <- ifelse(result$pvals<0.001,'***',ifelse(result$pvals<0.01,'**',ifelse(result$pvals<0.05,'*','')))
for (keyword in c('pval','FDR','lfc')){
    df1 = result[result$keyword==keyword,]
    df1$lab=NA
    if (keyword=='pval'){
        cutoffs=pvalcutoffs
        for (cutoff in cutoffs){
            df1[df1$cutoff==cutoff,'lab']=paste0('pvalue <= ',cutoff)
        }
    }else if (keyword=='FDR'){
        cutoffs=FDRcutoffs
        for (cutoff in cutoffs){
            df1[df1$cutoff==cutoff,'lab']=paste0('FDR <= ',cutoff)
        }
    }else{
        cutoffs=lfccutoffs
        for (cutoff in cutoffs){   
            df1[df1$cutoff==cutoff,'lab']=paste0('absLFC >= ', cutoff)
        }
    }
    p<-ggplot(data=df1,aes(y=celltype))+
        geom_point(aes(x=OR))+
        geom_segment(aes(x=lower,xend=upper,yend=celltype))+
        ylab('Cell Type') +
        xlab('Odds Ratio')+
        facet_wrap(~lab,nrow=length(cutoffs))+
        geom_vline(xintercept=1,linetype="dotted") +
        geom_text(aes(x=OR,label=plab),vjust=-0.5,fontface='bold')+
        scale_y_discrete(limits=rev)+
        theme_bw()+
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    pdf(paste0(output_path,'ATAC_OR_test_cutoff_',keyword,'.pdf'), useDingbats=FALSE,width=4,height=9)
    print(p)
    dev.off()
}
