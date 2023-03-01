library(stringr)
library(tidyverse)  
library(ggplot2)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/'
GSEA=read.table('/scratch/users/banwang3/celllines/scripts/RScripts/AGSEA_MWU_p-val_SignTest_Specific.csv',sep=',',header=T)
GSEA$celltype=str_split_fixed(GSEA$File, "_", 8)[,1]  
GSEA$db=str_split_fixed(GSEA$File, "_", 8)[,6]
GSEA$Up.in.Human=GSEA$Total-GSEA$Down.in.Human 
sigresult<-GSEA%>%group_by(celltype)%>% summarise(n = n())   
df=c()
for (ct in unique(GSEA$celltype)){
    CM=GSEA[GSEA$celltype==ct,]
    if (dim(CM)[1]>=5){
        CM=CM[0:5,]
    }
    tmp<-data.frame(term=factor(rep(c(1:nrow(CM)),2)),genes=c(-CM$Down.in.Human,CM$Up.in.Human),celltype=rep(CM$celltype,2),db=rep(CM$db,2),FDR=rep(CM$FDR.Sign.Test,2))
    tmp$Term=c(rep('',nrow(CM)),CM$Term) 
    df=rbind(df,tmp)
    #df$Term<-stringr::str_wrap(df$Term, 30)   
    #df$color <- ifelse(df$genes < 0, "Dodgerblue4", "Indianred3")
}
df$Term = str_split_fixed(df$Term, " - Homo", 2)[,1] 
df$Term = str_split_fixed(df$Term, ":R-HSA", 2)[,1] 
idx=which(apply(df, 1, function(x) any(grepl("GO:", x))))      
goterm=str_split_fixed(df[idx,'Term'],"~",2)[,2]
df[idx,'Term']=goterm
df$color <- ifelse(df$genes < 0, "dodgerblue2", "indianred1")  
celltypes=unique(GSEA$celltype)
celltypes=c('CM','MN')
for (ct in celltypes){
    pdf(paste0(output_path,ct,'_GSEA.pdf'),width=7,height=5,useDingbats=FALSE)
    df1=df[df$celltype==ct,]
    df1$dbterm=paste0(df1$db,': ',df1$Term) 
    df1$FDRlab <- ifelse(df1$FDR<0.01,'***',ifelse(df1$FDR<0.05,'**',ifelse(df1$FDR<0.25,'*','')))
    df1$dbterm[1:nrow(df1)/2]=''
    df1$FDRlab[1:nrow(df1)/2]='' 
    df1$Term = str_wrap(df1$Term,width=30) 
    df1$dbterm=firstup(df1$Term)
    print(ggplot(df1, aes(x = term, y = genes)) +
         geom_bar(stat = "identity", position = position_stack(),color = "white", fill = df1$color) +
        geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
        scale_x_discrete(limits=rev,breaks=c(1:(nrow(df1)/2)),labels=df1$dbterm[(nrow(df1)/2+1):nrow(df1)])+
        coord_flip()+
        theme_minimal()+
        ylab("Number of Genes")+
        theme(#axis.text.y = element_blank(),  # Remove Y-axis texts
            axis.text.y = element_text(size = 12),
            axis.title.y=element_blank(),
            axis.ticks.y = element_blank(), # Remove Y-axis ticks
            panel.grid.major.y = element_blank()) # Remove horizontal grid
        )
    dev.off()
}
idx = df[df$term %in% c(1:10),]
idx=which(apply(df, 1, function(x) any(grepl(c(1:10), x))))  

#Take the first 10 most significnat results.
nsig=10
subGSEA=GSEA[1:nsig,]
tmp<-data.frame(term=factor(rep(c(1:nrow(subGSEA)),2)),genes=c(-subGSEA$Down.in.Human,subGSEA$Up.in.Human),celltype=rep(subGSEA$celltype,2),db=rep(subGSEA$db,2),FDR=rep(subGSEA$FDR.Sign.Test,2))
tmp$Term=c(rep('',nrow(subGSEA)),subGSEA$Term) 
tmp$Term = str_split_fixed(tmp$Term, " - Homo", 2)[,1]
tmp$Term = str_split_fixed(tmp$Term, ":R-HSA", 2)[,1]
idx=which(apply(tmp, 1, function(x) any(grepl("GO:", x))))      
goterm=str_split_fixed(tmp[idx,'Term'],"~",2)[,2]
tmp[idx,'Term']=goterm
tmp$color <- ifelse(tmp$genes < 0, "dodgerblue2", "indianred1") 
tmp$dbterm=paste0(tmp$celltype,': ',tmp$Term,' (',tmp$db,')') 
tmp$FDRlab <- ifelse(tmp$FDR<0.01,'***',ifelse(tmp$FDR<0.05,'**',ifelse(tmp$FDR<0.25,'*','')))
tmp$dbterm[1:nrow(tmp)/2]=''
tmp$FDRlab[1:nrow(tmp)/2]='' 
tmp$dbterm=str_wrap(tmp$dbterm,width=40)
df1=tmp
pdf(paste0(output_path,'Top_',nsig,'_mostSig_GSEA.pdf'),useDingbats=FALSE)
print(ggplot(df1, aes(x = term, y = genes)) +
         geom_bar(stat = "identity", position = position_stack(),color = "white", fill = df1$color) +
         geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
         scale_y_continuous(limits = c(-max(abs(df1$genes)+1), max(abs(df1$genes)+1))) +
         geom_text(aes(label = FDRlab, hjust = -1, vjust = 0.5), size = 2.5) +
         scale_x_discrete(limits=rev,breaks=c(1:(nrow(df1)/2)),labels=df1$dbterm[(nrow(df1)/2+1):nrow(df1)])+
         coord_flip()+
         theme_minimal()+
         xlab("Term") +
         ylab("Number of Genes")+
         theme(#axis.text.y = element_blank(),  # Remove Y-axis texts
               axis.ticks.y = element_blank(), # Remove Y-axis ticks
               panel.grid.major.y = element_blank()) # Remove horizontal grid
         )
dev.off()
