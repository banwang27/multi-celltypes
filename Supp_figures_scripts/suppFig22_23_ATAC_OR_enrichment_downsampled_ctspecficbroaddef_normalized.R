library(stringr)
library(ggplot2)
library(dplyr)
library(epitools)
celltypes=c('CM','MN','HP','PP','SKM')
reffed='humreffed'
align2='human'
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Downsampled/BothAlleles/'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Downsampled/BothAlleles/'
allfiles = list.files(path=input_path,pattern='humreffed_all_peaks_down.txt')
df = read.csv(paste0(input_path, allfiles[1]),sep='\t',header=F)
colnames(df)=c('peak_id',allfiles[1])
for (filename in allfiles[2:length(allfiles)]){
    df1 = read.csv(paste0(input_path, filename),sep='\t',header=F)  
    colnames(df1)=c('peak_id',filename)
    df = merge(df,df1)
}
df = df[-seq(1,5),]
write.table(df,file=paste0(output_path,'AllSamples_BothAlleles_Counts_down.txt'),sep='\t',quote=F,row.names=F,col.names=T)
tmpdf = as.matrix(df[,-1])
tmpdf = tmpdf+1
CPM = sweep(as.matrix(tmpdf), 2, as.double(colSums(tmpdf)/1000000), `/`)
colnames(CPM)=paste0(str_split_fixed(colnames(CPM),'_BothAlleles',3)[,1],'_down_CPM')
CPM = data.frame(CPM)
CPM$peak_id=df$peak_id
write.table(CPM,file=paste0(output_path,'AllSamples_BothAlleles_CPM_down.txt'),sep='\t',quote=F,row.names=F,col.names=T)
idx = grep('PP',colnames(CPM))
PP = rowMeans(CPM[,idx])
idx = grep('MN',colnames(CPM))
MN = rowMeans(CPM[,idx])
idx = grep('CM',colnames(CPM))
CM = rowMeans(CPM[,idx])
CPM1 = cbind(CPM$peak_id,CM,CPM[,grep('HP',colnames(CPM))],MN,PP,CPM[,grep('SKM',colnames(CPM))])
CPM1 = data.frame(CPM1)
colnames(CPM1)=c('peak_id','CM','HP','MN','PP','SKM')
 write.table(CPM1,file=paste0(output_path,'Celltypes_BothAlleles_CPM_down.txt'),sep='\t',quote=F,row.names=F,col.names=T)
CPM = read.csv(paste0(output_path,'Celltypes_BothAlleles_CPM_down.txt'),sep='\t')
celltypes=c('CM','HP','MN','PP','SKM')
namestowrite=c()
final=c()
for (celltype in celltypes){
    for (ct1 in celltypes[-grep(celltype,celltypes)]){
        LFC = log2(CPM[,celltype]/CPM[,ct1])
        namestowrite = c(namestowrite,paste0(celltype,'_',ct1))
        final = cbind(final,LFC)
    }
}
final=data.frame(final)
colnames(final)=namestowrite
final$peak_id = CPM$peak_id
write.table(final,file=paste0(output_path,'Pairwise_LFC_BothAlleles_CPM_down.txt'),sep='\t',quote=F,row.names=F,col.names=T)
final = read.csv('/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Downsampled/BothAlleles/Pairwise_LFC_BothAlleles_CPM_down.txt',sep='\t')
celltypes=c('CM','HP','MN','PP','SKM')
result=c()
pvalcutoffs=c(0.1,0.05,0.01)
#lfccutoffs=c(0.5,1,1.5)
#input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Downsampled/'
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Downsampled/Downsampled_counts_persample/Downsampled_normalized/Downsampled_Normalized_binom/'
for (celltype in celltypes){
    dfopen = read.csv(paste0(input_path,celltype,'_Humreffed_Down_Filtered_Normed_Binom.txt'),sep='\t')
    #dfopen = dfopen[,c('Peak.Humreffed','l2fc_Humr','Human.binom_pval','Human.FDR')]    
    #colnames(dfopen)=c('peak_id','lfc','pval','FDR')
    dfopen = dfopen[,c('Peak','Binomial.p.value')]
    colnames(dfopen)=c('peak_id','pval')
    dfopen$celltype=celltype
    final1 = final[final$peak_id %in% dfopen$peak_id,] 
    for (i in seq(1,3)){
        idx1 = which(dfopen$pval<=pvalcutoffs[i])
        dCA_pval=dfopen[idx1,]
        tmpfinal=abs(final1[,grep(paste0(celltype,'_'),colnames(final1))])
        #idx_ctspecific = which(rowSums(tmpfinal > lfccutoffs[i]) ==(length(celltypes)-1))
        idx_ctspecific = which(rowSums(tmpfinal >=0.5) ==(length(celltypes)-1))
        dfopen_ctspecific = final1[idx_ctspecific,]
        OR<-c()
        lower<-c()
        upper<-c()
        pvals<-c()
        keywords<-c()
        cutoffs<-c()
        for (keyword in c('pval')){
            if (keyword=='pval'){
                dCA=dCA_pval
                cutoff=pvalcutoffs[i]
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
        print(keyword)
        print(cutoff)
        cutoffs<-c(cutoffs,cutoff)
        print(dat)
    }
    ORtest=data.frame(celltype=celltype,OR=OR,pvals=pvals,lower=lower,upper=upper,cutoff=cutoffs,keyword=keywords)
    result=rbind(result,ORtest)
}
}
result$plab <- ifelse(result$pvals<0.001,'***',ifelse(result$pvals<0.01,'**',ifelse(result$pvals<0.05,'*','')))
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/' 
for (keyword in c('pval')){
    df1 = result[result$keyword==keyword,]
    df1$lab=NA
    if (keyword=='pval'){
        cutoffs=pvalcutoffs
        for (cutoff in cutoffs){
            df1[df1$cutoff==cutoff,'lab']=paste0('pvalue <= ',cutoff)
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
    pdf(paste0(output_path,'ATAC_OR_test_Downsampled_cutoff_',keyword,'_pairwise_ctspecific_fix_absLFC0.5_normalized.pdf'), useDingbats=FALSE,width=4,height=9)
    print(p)
    dev.off()
}
    

