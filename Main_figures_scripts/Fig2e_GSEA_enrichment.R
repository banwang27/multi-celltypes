library(stringr)
library(dplyr)
library(ggplot2)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/'
GSEA=read.table('/scratch/users/banwang3/celllines/scripts/RScripts/AGSEA_MWU_p-val_SignTest_Specific_025_test.csv',sep=',',header=T)
GSEA$celltype=str_split_fixed(GSEA$File, "_", 8)[,1]  
GSEA$db=str_split_fixed(GSEA$File, "_", 8)[,6]
GSEA$Up.in.Human=GSEA$Total-GSEA$Down.in.Human 
GSEA=GSEA[order(GSEA$FDR.Sign.Test),]
sigresult<-GSEA%>%group_by(celltype)%>% summarise(n = n())    
df=c()
for (ct in unique(GSEA$celltype)){
    CM=GSEA[GSEA$celltype==ct,]
    if (dim(CM)[1]>=4){
        CM=CM[0:4,]
    }
    tmp<-data.frame(term=factor(rep(c(1:nrow(CM)),2)),genes=c(-CM$Down.in.Human,CM$Up.in.Human),celltype=rep(CM$celltype,2),db=rep(CM$db,2),FDR=rep(CM$FDR.Sign.Test,2))
    tmp$Term=c(rep('',nrow(CM)),CM$Term) 
    df=rbind(df,tmp)
}
df$Term = str_split_fixed(df$Term, " - Homo", 2)[,1] 
df$Term = str_split_fixed(df$Term, ":R-HSA", 2)[,1] 
idx=which(apply(df, 1, function(x) any(grepl("GO:", x))))      
goterm=str_split_fixed(df[idx,'Term'],"~",2)[,2]
df[idx,'Term']=goterm 
df$color <- ifelse(df$genes < 0, "#1E90FF","#FF6A6A")
celltypes=unique(GSEA$celltype)
celltypes=c('CM','MN')
for (ct in celltypes){
    pdf(paste0(output_path,ct,'_GSEA.pdf'),width=5.5,height=4,useDingbats=FALSE)
    df1=df[df$celltype==ct,]
    maxvalue = max(abs(df1$genes))
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
        ylim(-maxvalue,maxvalue)+
        xlab("")+
        theme(#axis.text.y = element_blank(),  # Remove Y-axis texts
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 15),
            axis.ticks.x=element_line(colour="black"),
            axis.line.x = element_line(colour="black"),
            axis.title.x=element_text(size=15),
            axis.title.y=element_blank(),
            axis.ticks.y = element_blank(), # Remove Y-axis ticks
            panel.grid.major.y = element_blank() # Remove horizontal grid
        ) 
        )
    dev.off()
}
