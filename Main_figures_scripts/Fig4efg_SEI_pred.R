library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
getCorString <- function(dataframe, x, y, category_col = NULL, labelRound = 2, method = "pearson") {
    # Construct category list of factors
    if(!is.null(category_col)) {
        catList <- list()
        for(i in category_col){
            catList[[i]] <- dataframe[,i]
           }
    } else {
         catList <- rep(T, nrow(dataframe))
    }
    # Loop by categorie
    results <- by(dataframe, catList, function(x, xVal, yVal, category_col) { 
                  # Get correlations
                  corR <- cor.test(x[,xVal], x[,yVal], method = method)
                  RMSE=sqrt(mean((x[,yVal]-x[,xVal])**2))
                  #corText <- paste0("r = ", signif(corR$estimate,labelRound) , "\np = ", signif(corR$p.value, labelRound),"\nslope = ", signif(coef(lmtest)[2],labelRound))
                  corText <- paste0("RMSE = ",signif(RMSE,labelRound))
                  result <- data.frame(r = corR$estimate, p =  corR$p.value, text = corText, stringsAsFactors = F)
                  # Append category names
                  if(!is.null(category_col)) {
                      for(i in category_col){
                       result <- cbind (result, x[1,i], stringsAsFactors = F)
                      }
                      colnames(result)[ (ncol(result) - length(category_col) + 1) : ncol(result) ] <- category_col
                   }
                   return(result)
                }, xVal = x, yVal = y, category_col = category_col)
     results <- do.call(rbind, results)
     return(results)
}
plotSEIrawoutput<-function(df1,xaxis,yaxis,idcandidate=c('Neuron.H3K4me3.ID.64664','iPSC.H3K4me3.ID.75085'),xlim1=0.2,xlim2=1,ylim1=0.2,ylim2=1){
     data = data.frame(x=df1[,xaxis],y=df1[,yaxis],id=df1[,'id'])  
     pointlabel=data[data$id %in% idcandidate,]
     pointlabel$id=str_split_fixed(pointlabel$id,'.ID',2)[,1]
     pointlabel$id1=str_split_fixed(pointlabel$id,'\\.',2)[,1]
     pointlabel$id2=str_split_fixed(pointlabel$id,'\\.',2)[,2]
     pointlabel$id=paste0(pointlabel$id1,'\n',pointlabel$id2)
     texts=getCorString(data,'x','y')  
     texts$x=-Inf   
     texts$y=Inf
     p<- ggplot(data,aes(x=x,y=y))+     
         geom_hline(yintercept=0, color = "grey")+
         geom_vline(xintercept=0, color = "grey")+
         #geom_smooth(method=lm,se=FALSE, fullrange=TRUE,linewidth=0.5,alpha=0.2,linetype = "longdash")+ 
         geom_point(size=0.2,alpha=0.3,show.legend = FALSE)+
         geom_text(data=texts,aes(label=text),hjust=0,vjust=1,size=3)+
         xlab(xaxis)+
         ylab(yaxis)+
         xlim(xlim1,xlim2)+
         ylim(ylim1,ylim2)+
        theme_bw()+ 
         theme(axis.title = element_text(size = 8),legend.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    p<-p+
        geom_label_repel(data=pointlabel,aes(x=x,y=y,label=id),box.padding   = 0.35, segment.color = 'black',size=3,min.segment.length = 0)+
        geom_point(data=pointlabel,aes(x=,y=y),color="red",size=1)
     return(p)
}
plotSEIhistogram<-function(df1,xaxis,pointlabel,xlim1=0.2,xlim2=1,ylim1=0,ylim2=1){
    data = data.frame(x=df1[,xaxis])
    p<-ggplot(data,aes(x=x))+
        geom_hline(yintercept=0, color = "black")+
        geom_density(alpha=0.2,fill="white",colour="black")+
        xlab(xaxis)+
        xlim(xlim1,xlim2)+   
        ylim(ylim1,ylim2)+
        ylab('')+
        theme_bw()+
        theme(axis.title = element_text(size = 8),legend.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    
    toplot=pointlabel
    toplot$y=0
    p<-p+
        geom_label_repel(data=toplot,aes(x=score,y=y,label=id),box.padding   = 0.35, segment.color = 'black',size=3,min.segment.length = 0)+
        geom_point(data=toplot,aes(x=score,y=y),color="red",size=1)
   return(p)
}
makePointlabel<-function(df,xaxis,idcandidate=c('Neuron.H3K4me3.ID.64664','iPSC.H3K4me3.ID.75085')){
     df4 = df[df$id %in% idcandidate,]
     pointlabel=data.frame(score=df4[,xaxis],id=df4[,'id'])
     pointlabel$id=str_split_fixed(pointlabel$id,'.ID',2)[,1]
     pointlabel$id1=str_split_fixed(pointlabel$id,'\\.',2)[,1]
     pointlabel$id2=str_split_fixed(pointlabel$id,'\\.',2)[,2]
     pointlabel$id=paste0(pointlabel$id1,'\n',pointlabel$id2)
     return(pointlabel[,c('id','score')])
}



input_path='/scratch/users/banwang3/references/outputs/FABP7_centeratpromoter/sei/'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/'

df = read.csv('FABP7_centeratpromoter_summary_SEI.txt',sep='\t')
idx = which(df[,paste0('humanfasta')]>=0.5| df[,paste0('chimpfasta')]>=0.5)
df1 = df[idx,]
print(nrow(df1))

p1<-plotSEIhistogram(df1,'humanfasta',makePointlabel(df1,'humanfasta'),ylim1=0,ylim2=2.5)
p2<-plotSEIhistogram(df1,'chimpfasta',makePointlabel(df1,'chimpfasta'),ylim1=0,ylim2=2.5)
p3<-plotSEIhistogram(df1,'humanfasta_withindel',makePointlabel(df1,'humanfasta_withindel'),ylim1=0,ylim2=2.5)
p4<-plotSEIhistogram(df1,'humanfasta_withSNP26',makePointlabel(df1,'humanfasta_withSNP26'),ylim1=0,ylim2=2.5)
p<-(p1|p4)/(p2|p3)
pdf(paste0(output_path,'FABP7_SEI_hist.pdf'),useDingbats=FALSE,width=6,height=6)
print(p)
dev.off()

ids = c('hum_indelhum_fastalfc','hum_humdel1G_fastalfc','hum_humadd3bp_fastalfc','humanbased_humchplfc_withSNP26')
for (id in ids){
    p5<-plotSEIrawoutput(df1,id,'hum_chp_fastalfc',ylim1=-0.54,ylim2=1.42,xlim1=-0.54,xlim2=1.42) 
    pdf(paste0(output_path,'FABP7_SEI_lfc_',id,'.pdf'),width=3,height=3)
    print(p5)
    dev.off()
}
p5<-plotSEIrawoutput(df1,'hum_chp_fastalfc','hum_indelhum_fastalfc',xlim1=-0.4,xlim2=1.42,ylim1=-0.4,ylim2=1.42)
p6<-plotSEIrawoutput(df1,'hum_chp_fastalfc','humanbased_humchplfc_withSNP1',xlim1=-0.4,xlim2=1.42,ylim1=-0.4,ylim2=1.42)  
p<-p6/p5
pdf(paste0(output_path,'FABP7_SEI_lfc.pdf'),useDingbats=FALSE,width=3,height=6)
print(p)
dev.off()
#Human/Chimp based single/All SNP flipping result
constructBoxplot<-function(df1,species,includeSNPall=T){
    idx = grep(paste0(species,'based_humchplfc_withSNP'),colnames(df1))
    SNPname=colnames(df1)[idx[length(idx)]]
    nSNP=as.numeric(str_split_fixed(SNPname,'withSNP',2)[,2])
    df2=c()
    for (i in seq(1,nSNP)){  
        sei = data.frame(lfc = df1[,paste0(species,'based_humchplfc_withSNP',i)],SNP = paste0('SNP',i),type='SNP effect') 
        df2=rbind(df2,sei)
    }
    if (species=='human'){
        allSNPflipid='hum_hum2chp_fastalfc'
    }else{
        allSNPflipid='chp2hum_chp_fastalfc'
    }
    sei=data.frame(lfc = df1[,'hum_chp_fastalfc'],SNP='Chimp',type='Net effect')
    df2 = rbind(df2,sei)
    sei=data.frame(lfc=df1[,'hum_humdel1G_fastalfc'],SNP='Deletion',type='Individual candidate causal variant effect')
    df2 = rbind(df2,sei)
    sei=data.frame(lfc=df1[,'hum_humadd3bp_fastalfc'],SNP='Insertion',type='Individual candidate causal variant effect')
    df2 = rbind(df2,sei)
    sei=data.frame(lfc=df1[,'hum_indelhum_fastalfc'],SNP='Insertion\n+Deletion',type='Combined candidate causal variant effect')
    df2 = rbind(df2,sei)
    sei=data.frame(lfc = df1[,allSNPflipid],SNP='SNPall',type='Net effect')
    df2 = rbind(df2,sei)
    df2$SNP = factor(df2$SNP,levels=c('Chimp','Insertion\n+Deletion','Deletion','Insertion','SNPall',paste0('SNP',seq(1,nSNP))))
    df2$type = factor(df2$type,levels=c('Net effect','Combined candidate causal variant effect','Individual candidate causal variant effect','SNP effect'))
    if (!includeSNPall){
        df2 = df2[df2$SNP!='SNPall',]
    }
    p1<-ggplot(df2,aes(x = SNP,y=lfc,fill=type))+
        geom_hline(yintercept=0, color = "grey")+
        geom_boxplot(outlier.size = 0.15,show.legend=T)+
        scale_fill_manual(values=c("#555599","#DD4444","#EB9090","#00BFC4"))+
        ylab('Log2(Human/non-Human)')+
        xlab('')+
        theme_bw()+
        theme(axis.text=element_text(size=14),axis.title=element_text(size=18))+
        theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    return(p1)
}
species='human'
p<-constructBoxplot(df1,species,includeSNPall=F)   
pdf(paste0(output_path,'FABP7_',species,'based_SNPeffects.pdf'),useDingbats=FALSE,width=15,height=5.5)
print(p)
dev.off()
for (species in c('human','chimp')){
    p<-constructBoxplot(df1,species,includeSNPall=F)
    pdf(paste0(output_path,'FABP7_',species,'based_SNPeffects.pdf'),useDingbats=FALSE,width=15,height=4.5)
    print(p)
    dev.off()
}
