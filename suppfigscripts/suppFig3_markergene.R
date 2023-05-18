library(ggplot2)
library(stringr)
backup_path='/scratch/users/banwang3/celllines/data_backup'
df = read.csv(paste0(backup_path,'/scratch/users/banwang3/celllines/analysis/ASE_norm/human/filtered/TPM.txt'),sep='\t')
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/'
idx =grep('total',colnames(df))
idx1=grep('gene',colnames(df))
df = df[,c(idx1,idx)]
celltypes = c('iPSC','CM','HP','MN','SKM','PP','RPE')
targetgenedf=c()
for (celltype in celltypes){
    if (celltype=='CM'){  
        subdf<-data.frame(celltype=celltype,target_genes=c('TNNT2', 'NKX2-5', 'MYBPC3') )
    }else if(celltype=='HP'){  
        subdf<-data.frame(celltype=celltype,target_genes=c('ALB','TTR','GPC3','HNF1B'))
    }else if(celltype=='MN'){
        subdf<-data.frame(celltype=celltype,target_genes=c('OLIG2', 'ELAVL2', 'ISL1'))
    }else if (celltype=='SKM'){
        subdf<-data.frame(celltype=celltype,target_genes=c('MYLPF', 'MYOG','MYL1'))
    }else if (celltype=='PP'){
        subdf<-data.frame(celltype=celltype,target_genes=c('PDX1', 'NKX6-1','GP2','NEUROG3'))
    }else if (celltype == 'RPE'){
        subdf<-data.frame(celltype=celltype,target_genes=c('PAX6','TYRP1','MITF'))
    }else if (celltype=='iPSC'){
        subdf<-data.frame(celltype=celltype,target_genes=c('NANOG','MKI67'))
    }
    targetgenedf = rbind(targetgenedf,subdf)
}
for (celltype in celltypes){
    idx = grep(celltype,colnames(df))
    tmp = cbind(tmp,rowMeans(df[,idx]))
}
colnames(tmp)=celltypes
tmp = data.frame(tmp)
tmp$gene=df$gene
df1 = df[df$gene %in% targetgenedf$target_genes,]
genenames=df1$gene
target_genes=genenames
idx = grep('gene',colnames(df1))
sampleTPM = c()
for (gene in genenames){
    tmp1 = t(df1[df1$gene==gene,-idx]) 
    tmp1=as.data.frame(tmp1)
    colnames(tmp1)='TPM' 
    tmp1$sampleid=row.names(tmp1)
    tmp1$gene=gene
    sampleTPM=rbind(sampleTPM,tmp1)
}
sampleTPM$celltype=str_split_fixed(sampleTPM$sampleid,'_',4)[,2]
sampleTPM$hyb=str_split_fixed(sampleTPM$sampleid,'_',4)[,3]

sampleTPM$celltype=factor(sampleTPM$celltype)
sampleTPM$hyb=ifelse(sampleTPM$hyb=='Hyb1','Hybrid1','Hybrid2')
sampleTPM<-merge(sampleTPM,targetgenedf,by.x='gene',by.y='target_genes',all.x=TRUE,suffixes=c('','_marker'))
pboxplot<-function(sampleTPM,celltype_marker){
    sampleTPM<-sampleTPM[sampleTPM$celltype_marker==celltype_marker,]
    p<-ggplot(data=sampleTPM,aes(x=celltype,y=TPM))+
        geom_boxplot(outlier.size=0,fill="white",colour="black",width=0.5)+
        geom_point(data=sampleTPM,aes(x=celltype,y=TPM,fill=hyb), size=1.5, shape=21, colour="grey20",position=position_jitter(width=0.1, height=0.1))+
        scale_fill_manual(values=c("#6959cd","#3cb371")) +
        
        facet_wrap(~gene,ncol=length(unique(sampleTPM$gene)),scales = "free")+ 
        #facet_wrap(~gene,scales="free")+
        xlab('')+
        ylab('Expression level (TPM)')+
        ggtitle(paste0('Marker genes for ',celltype_marker))+
        theme_bw()+    
        guides(fill=guide_legend(title="Hybrid line"))+
        theme(strip.text = element_text(face = "italic"))+
        theme(legend.position="top",panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    return(p)
}
printp<-function(p,celltype,width){
    if (celltype!='iPSC'){
        p<-p+theme(legend.position="none")
    }
    pdf(paste0(output_path,celltype,'_TPM_markergene_final.pdf'),useDingbats=FALSE,width=width,height=3)
    print(p)
    dev.off()
}

for (celltype_marker in celltypes){
    p1<-pboxplot(sampleTPM,celltype_marker)
    
    widths=length(unique(sampleTPM[sampleTPM$celltype_marker==celltype_marker,'gene']))*2.5 
    printp(p1,celltype_marker,widths)
}

#Main figure
final_genes=c('NKX2-5','HNF1B','OLIG2','GP2','PAX6','MYLPF')
sampleTPM1 = sampleTPM[sampleTPM$gene %in% final_genes,]
sampleTPM1$gene=factor(sampleTPM1$gene, levels=final_genes)
pdf(paste0(output_path,'marker_gene_TPM.pdf'),useDingbats=FALSE,width=7,height=5)
p<-pboxplot(sampleTPM1)
print(p)
dev.off()
