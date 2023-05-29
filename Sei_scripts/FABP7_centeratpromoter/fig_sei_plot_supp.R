library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
getCorString <- function(dataframe, x, y, category_col = NULL, labelRound = 2, method = "pearson") {
    # Calculates pearson correlation coefficient and pvalue between two columns in a dataframe
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
                  lmtest<-lm(x[,yVal]~x[,xVal])
                  #corText <- paste0("r = ", signif(corR$estimate,labelRound) , "\np = ", signif(corR$p.value, labelRound),"\nslope = ", signif(coef(lmtest)[2],labelRound))
                  corText <- paste0("R\u00B2 = ",signif(summary(lmtest)$r.squared,labelRound))
                  #corText <- paste0("MSE = ",signif(mean((x[,xVal] - x[,yVal])^2),labelRound),"\nR^2 = ",signif(summary(lmtest)$r.squared,labelRound))
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
    p<-p+geom_point(data=pointlabel,aes(x=,y=y),color="red",size=0.2)+
        geom_label_repel(data=pointlabel,aes(x=x,y=y,label=id),box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50',size=3) 
        
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
    #pp<-ggplot_build(p)
    #newdata=data.frame(x = pp$data[[1]]$x,y=pp$data[[1]]$density)
    #toplot=c()
    #for (i in seq(1,nrow(pointlabel))){
    #    toplot=rbind(toplot,newdata[which.min(abs(newdata$x-pointlabel$score[i])),])
    #}
    #toplot$id=pointlabel$id
    toplot=pointlabel
    toplot$y=0
    p<-p+
        geom_point(data=toplot,aes(x=score,y=y),color="red",size=1)+
        geom_label_repel(data=toplot,aes(x=score,y=y,label=id),box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50',size=3) 
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
output_path='/scratch/users/banwang3/references/scripts/seifigures/'


df = read.csv('FABP7_centeratpromoter_summary_SEI.txt',sep='\t')
idx = which(df[,paste0('humanfasta')]>=0.5| df[,paste0('chimpfasta')]>=0.5)
df1 = df[idx,]
print(nrow(df1))

humanfasta_del1G humanfasta_add3bp
p1<-plotSEIhistogram(df1,'humanfasta_del1G',makePointlabel(df1,'humanfasta_del1G'),ylim1=0,ylim2=2.5)
p2<-plotSEIhistogram(df1,'humanfasta_add3bp',makePointlabel(df1,'humanfasta_add3bp'),ylim1=0,ylim2=2.5)
#p3<-plotSEIhistogram(df1,'humanfasta_withindel',makePointlabel(df1,'humanfasta_withindel'),ylim1=0,ylim2=2.5)
#p4<-plotSEIhistogram(df1,'humanfasta_withSNP1',makePointlabel(df1,'humanfasta_withSNP1'),ylim1=0,ylim2=2.5)
#p<-(p1|p4)/(p2|p3)
#pdf(paste0(output_path,'FABP7_SEI_hist.pdf'),useDingbats=FALSE,width=6,height=6)
#print(p)
#dev.off()

p5<-plotSEIrawoutput(df1,'hum_chp_fastalfc','hum_humdel1G_fastalfc',xlim1=-0.4,xlim2=1.42,ylim1=-0.4,ylim2=1.42)
p6<-plotSEIrawoutput(df1,'hum_chp_fastalfc','hum_humadd3bp_fastalfc',xlim1=-0.4,xlim2=1.42,ylim1=-0.4,ylim2=1.42)  
p<-(p1|p5)/(p2|p6)
pdf(paste0(output_path,'FABP7_SEI_supp_del1G_add3bp.pdf'),useDingbats=FALSE,width=6,height=6)
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
    sei=data.frame(lfc = df1[,allSNPflipid],SNP='SNPall',type='Net effect')
    df2 = rbind(df2,sei)
    df2$SNP = factor(df2$SNP,levels=c('Chimp','SNPall',paste0('SNP',seq(1,nSNP))))
    df2$type = factor(df2$type,levels=c('Net effect','SNP effect'))
    if (!includeSNPall){
        df2 = df2[df2$SNP!='SNPall',]
    }
    p1<-ggplot(df2,aes(x = SNP,y=lfc,fill=type))+
        geom_hline(yintercept=0, color = "grey")+
        geom_boxplot(outlier.size = 0.15,show.legend=T)+
        scale_fill_manual(values=c("#555599","#00BFC4"))+
        ylab('Log2(Human/non-Human)')+
        xlab('')+
        theme_bw()+
        theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    return(p1)
}
for (species in c('human','chimp')){
    p<-constructBoxplot(df1,species,includeSNPall=F)
    pdf(paste0(output_path,'FABP7_',species,'based_SNPeffects.pdf'),useDingbats=FALSE,width=15,height=4.5)
    print(p)
    dev.off()
}
