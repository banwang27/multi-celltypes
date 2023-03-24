df = read.csv('/scratch/users/banwang3/celllines/analysis/ASE_norm/human/filtered/ASE_value_fromrawcounts.txt',sep='\t')
rawcounts = read.csv('/scratch/users/banwang3/celllines/analysis/ASE_norm/human/filtered/rawcounts.txt',sep='\t')
idxa=grep('total',colnames(rawcounts))
rawcounts=rawcounts[,-idxa]
genenamesall=df$gene
idx1 = grep('gene',colnames(df))
df = df[,-idx1]
for (includeSKM in c(T,F)){
  if (includeSKM){
    celltypes=c('CM','MN','HP','PP','RPE','SKM')
    SKM='includeSKM'
  }else{
    celltypes=c('CM','MN','HP','PP','RPE')
    idx1 = grep('SKM',colnames(df))
    df = df[,-idx1]
    SKM='excludeSKM'
  }
  for (celltype in celltypes){
    idx=which(rowSums(abs(df))>0,)
    revisedf=df[idx,]
    rawcounts1=rawcounts[idx,grep(celltype,colnames(rawcounts))]
    humallele = rowMeans(rawcounts1[,grep('hum',colnames(rawcounts1))])
    chpallele = rowMeans(rawcounts1[,grep('chp',colnames(rawcounts1))])
    ASEvalue=revisedf[,celltype]
    genenames=genenamesall[idx]
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
    sumofallgenes = colSums(expr)
    a=sum(sumofallgenes)/sumofallgenes[celltype]
    b=expr[,celltype]/rowSums(expr)
    PEM=a*b
    PEM = PEM/max(PEM)
    PEM=data.frame(gene=genenames,dEE=PEM,ASE=ASEvalue,humallele=humallele,chpallele=chpallele)
    write.table(PEM,file=paste0('TSmetrics/EE_fromrawcounts/dEE_zeroing_for',celltype,'_',SKM,'.txt'),quote=F,sep='\t',row.names=F,col.names=T)
  }
}
