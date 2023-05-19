library(stringr)
input_path='/scratch/users/banwang3/celllines/Enrichment_test/ATAC_fortau/Unfiltered/'
allfiles = list.files('/scratch/users/banwang3/celllines/Enrichment_test/ATAC_fortau/Unfiltered/',pattern='.txt')

includeSKM=T
if (includeSKM){
    celltypes=c('CM','MN','HP','PP','SKM')
    SKM='includeSKM'
}else{
    celltypes=c('CM','MN','HP','PP')
    SKM='excludeSKM'
}
df1=c()
df2=c()
allelecounts=c()
counts_humr=c()
allelecounts_chp=c()
counts_chpr=c()
for (celltype in celltypes){
    df = read.csv(paste0(input_path,allfiles[grep(celltype,allfiles)]),sep='\t')
    df = df[order(df$Peak.Humreffed),]
    dftmp=df[,grep('humreffed',colnames(df))]
    allelecounts=c(allelecounts,sum(colSums(dftmp))/(ncol(dftmp)/2))
    counts_humr=cbind(counts_humr,rowSums(dftmp)/(ncol(dftmp)/2))
    
    dftmp=df[,grep('chpreffed',colnames(df))]   
    allelecounts_chp=c(allelecounts_chp,sum(colSums(dftmp))/(ncol(dftmp)/2))
    counts_chpr = cbind(counts_chpr,rowSums(dftmp)/(ncol(dftmp)/2)) 
    df1=cbind(df1,df$l2fc_Humr)
    df2=cbind(df2,df$l2fc_Chpr)
}
allelecounts=data.frame(celltype=celltypes,allelecounts_humr=allelecounts,allelecounts_chpr=allelecounts_chp)
allelecounts$factor=allelecounts$allelecounts_humr/min(allelecounts$allelecounts_humr)
colnames(df1)=paste0(celltypes,'_lfc')
colnames(df2)=paste0(celltypes,'_lfc')
colnames(counts_humr)=paste0(celltypes,'_allelecounts')
colnames(counts_chpr)=paste0(celltypes,'_allelecounts')
humr = cbind(df1,counts_humr)
humr = data.frame(humr)
humr$peak_id=df$Peak.Humreffed
chpr = cbind(df2,counts_chpr)
chpr = data.frame(chpr)
chpr$peak_id = df$Peak.Chpreffed
counts=humr[,grep('allelecounts',colnames(humr))]
idx = which(rowSums(counts>10)==ncol(counts))
humr=humr[idx,]
counts=humr[,grep('allelecounts',colnames(humr))]
allelecounts=allelecounts[,c('celltype','factor')]
idxtocount=c()

counts=humr[,grep('allelecounts',colnames(humr))]


computeready<-function(revisedf,celltype){
    targetgenelist = revisedf[,celltype]
    idx = grep(celltype,colnames(revisedf)) 
    testdf = revisedf[,-idx]
    tmpdf = targetgenelist*testdf
    for (i in seq(1,ncol(tmpdf))){
      idx = which(tmpdf[,i]<0)
      celltypename = colnames(tmpdf)[i] 
      revisedf[idx,celltypename]=0
    }
    expr = abs(revisedf)
    expr = data.frame(expr)
    return(expr)
}

            
computedEE<-function(expr,celltype){
        sumofallgenes = colSums(expr)
        a = sum(sumofallgenes)/sumofallgenes[celltype]
        b=expr[,celltype]/rowSums(expr)
        EE=a*b
        EE=EE/max(EE)
        return(EE)
}

colnames(df1)=celltypes
colnames(df2)=celltypes
for (celltype in celltypes){
    expr = computeready(df1,celltype)
    expr_chp=computeready(df2,celltype)
    
    EE = computedEE(expr,celltype)
    EE_chp=computedEE(expr_chp,celltype)
    
    dEE=data.frame(Peak.Humreffed=df$Peak.Humreffed,dEE.Humreffed=EE,Peak.Chpreffed=df$Peak.Chpreffed,dEE.Chpreffed=EE_chp)
    write.table(dEE,file=paste0('dEE_zeroing_for',celltype,'_',SKM,'.txt'),quote=F,sep='\t',row.names=F,col.names=T) 
}


