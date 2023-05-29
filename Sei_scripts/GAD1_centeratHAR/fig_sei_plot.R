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
                  RMSE=sqrt(mean((x[,yVal]-x[,xVal])**2))
                  #corText <- paste0("r = ", signif(corR$estimate,labelRound) , "\np = ", signif(corR$p.value, labelRound),"\nslope = ", signif(coef(lmtest)[2],labelRound))
                  #corText <- paste0("MSE = ",signif(mean((x[,xVal] - x[,yVal])^2),labelRound))
                  #corText <- paste0("R\u00B2 = ",signif(summary(lmtest)$r.squared,labelRound))
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
plotSEIrawoutput<-function(df1,xaxis,yaxis,idcandidate=c('iPSC.H3K4me3.ID.75085','K562_Erythroblast_Bone_Marrow.ZNF639.ID.64649'),xlim1=0.2,xlim2=1,ylim1=0.2,ylim2=1){
     data = data.frame(x=df1[,xaxis],y=df1[,yaxis],id=df1[,'id'])  
     pointlabel=data[data$id %in% idcandidate,]
     pointlabel$id = ifelse(grepl('iPSC',pointlabel$id),'iPSC\nH3K4me3','K562\nZNF639')
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
        #geom_point(data=pointlabel,aes(x=,y=y),color="red",size=0.2)+
        #geom_label_repel(data=pointlabel,aes(x=x,y=y,label=id),box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50',size=3)
        geom_label_repel(data=pointlabel,aes(x=x,y=y,label=id),box.padding   = 0.35, segment.color = 'black',size=3,min.segment.length = 0)+ 
        geom_point(data=pointlabel,aes(x=,y=y),color="red",size=1)  
        
     return(p)
}
plotSEIhistogram<-function(df1,xaxis,pointlabel,xlim1=0,xlim2=1,ylim1=0,ylim2=2.5){
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
    toplot=pointlabel
    toplot$y=0
    #toplot$id=pointlabel$id
    p<-p+
        #geom_point(data=toplot,aes(x=score,y=y),color="red",size=1)+
        #geom_label_repel(data=toplot,aes(x=score,y=y,label=id),box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50',size=3)
        geom_label_repel(data=toplot,aes(x=score,y=y,label=id),box.padding   = 0.35, segment.color = 'black',size=3,min.segment.length = 0)+
        geom_point(data=toplot,aes(x=score,y=y),color="red",size=1)
   return(p)
}
makePointlabel<-function(df,xaxis,idcandidate=c('iPSC.H3K4me3.ID.75085','K562_Erythroblast_Bone_Marrow.ZNF639.ID.64649')){     
     df4 = df[df$id %in% idcandidate,]
     pointlabel=data.frame(score=df4[,xaxis],id=df4[,'id'])
     pointlabel$id = ifelse(grepl('iPSC',pointlabel$id),'iPSC\nH3K4me3','K562\nZNF639')
     return(pointlabel[,c('id','score')])
}



input_path='/scratch/users/banwang3/references/outputs/GAD1_centeratHAR/sei/'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/'

df = read.csv('GAD1_centeratHAR_summary_SEI.txt',sep='\t')
idx = which(df[,paste0('humanfasta')]>=0.5| df[,paste0('chimpfasta')]>=0.5)
df1 = df[idx,]
print(nrow(df1))

df1 = df1[order(df1$hum_chp_fastalfc),] 
ids = c("hum_chp_fastalfc","humanbased_humchplfc_withSNP1","humanbased_humchplfc_withSNP21")
print(df1[1:10,c('id',ids)])
tmp=df1[df1$hum_chp_fastalfc<0.01&df1$hum_chp_fastalfc>-0.01,] 
print(tmp[,c('id',ids)])   

p1<-plotSEIhistogram(df1,'humanfasta',makePointlabel(df1,'humanfasta'),ylim1=0,ylim2=2.5)
p2<-plotSEIhistogram(df1,'chimpfasta',makePointlabel(df1,'chimpfasta'),ylim1=0,ylim2=2.5)
p3<-plotSEIhistogram(df1,'humanfasta_withSNP1',makePointlabel(df1,'humanfasta_withSNP1'),ylim1=0,ylim2=2.5)
p4<-plotSEIhistogram(df1,'humanfasta_withSNP21',makePointlabel(df1,'humanfasta_withSNP21'),ylim1=0,ylim2=2.5)
p<-(p1|p3)/(p2|p4)
pdf(paste0(output_path,'GAD1_SEI_hist.pdf'),useDingbats=FALSE,width=6,height=6)
print(p)
dev.off()

ids=c('humanbased_humchplfc_withSNP1','humanbased_humchplfc_withSNP20','humanbased_humchplfc_withSNP21','humanbased_humchplfc_withSNP34')
for (id in ids){
    p5<-plotSEIrawoutput(df1,id,'hum_chp_fastalfc',xlim1=-2.7,xlim2=0.96,ylim1=-2.7,ylim2=0.96)
    pdf(paste0(output_path,'GAD1_SEI_lfc_',id,'.pdf'),width=3,height=3)
    print(p5)
    dev.off()
}
p6<-plotSEIrawoutput(df1,'hum_chp_fastalfc','humanbased_humchplfc_withSNP1',xlim1=-2.7,xlim2=0.7,ylim1=-2.7,ylim2=0.7) 
p<-p6/p5
pdf(paste0(output_path,'GAD1_SEI_lfc.pdf'),useDingbats=FALSE,width=3,height=6)
print(p)
dev.off()
#Human/Chimp based single/All SNP flipping result
constructBoxplot<-function(df1,species,includeSNPall=T){
    idx = grep(paste0(species,'based_humchplfc_withSNP'),colnames(df1))
    SNPname=colnames(df1)[idx[length(idx)]]
    nSNP=as.numeric(str_split_fixed(SNPname,'withSNP',2)[,2])
    df2=c()
    for (i in seq(1,nSNP)){  
        if (i==21){
            sei = data.frame(lfc = df1[,paste0(species,'based_humchplfc_withSNP',i)],SNP = paste0('SNV',i),type='Candidate causal variant effect')
        }else{
            sei = data.frame(lfc = df1[,paste0(species,'based_humchplfc_withSNP',i)],SNP = paste0('SNV',i),type='SNV effect') 
        }
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
    df2$SNP = factor(df2$SNP,levels=c('Chimp','SNPall',paste0('SNV',seq(1,nSNP))))
    #df2$type = factor(df2$type,levels=c('Net effect','SNP effect'))
    df2$type = factor(df2$type,levels= c('Net effect','Candidate causal variant effect','SNV effect'))
    if (!includeSNPall){
         df2 = df2[df2$SNP!='SNPall',]
    }
    p1<-ggplot(df2,aes(x = SNP,y=lfc,fill=type))+
        geom_hline(yintercept=0, color = "grey")+
        geom_boxplot(outlier.size = 0.15,show.legend=T)+
        scale_fill_manual(values=c("#555599","#DD4444","#00BFC4"))+
        ylab('Log2(Variant effects)')+
        xlab('')+
        theme_bw()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20))+
        theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    return(p1)
}
species='human'
p<-constructBoxplot(df1,species,includeSNPall=F) 
pdf(paste0(output_path,'GAD1_',species,'based_SNPeffects.pdf'),useDingbats=FALSE,width=15,height=5.5)
print(p)
dev.off()
for (species in c('human','chimp')){
    p<-constructBoxplot(df1,species,includeSNPall=F)
    pdf(paste0(output_path,'GAD1_',species,'based_SNPeffects.pdf'),useDingbats=FALSE,width=15,height=4.5)
    print(p)
    dev.off()
}
