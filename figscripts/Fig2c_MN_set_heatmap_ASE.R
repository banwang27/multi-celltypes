library(ggplot2)  
backup_path='/scratch/users/banwang3/celllines/data_backup'
align2='human'
TGF = c('NGF','BDNF','NTF3','NTF4','NTRK1','NTRK2','NTRK3','NGFR')
length(TGF)
TGF = data.frame('gene'=TGF)
final=c()
for (celltype in c('CM','MN','PP','HP','SKM','RPE')){
    deseq2 = read.csv(paste0(backup_path,'/scratch/users/banwang3/celllines/analysis/DESeq2_norm/human/',celltype,'_aligned2human_qvalue_DESeq2.txt'),sep='\t')
    deseq2m = read.csv(paste0(backup_path,'/scratch/users/banwang3/celllines/analysis/DESeq2_norm/merged/',celltype,'_merged_qvalue_DESeq2.txt'),sep='\t')
    deseq2$flag = deseq2$gene %in% deseq2m$gene
    deseq2=deseq2[,c('gene','log2FoldChange','qval_mine','flag')]
    colnames(deseq2)=c('gene',paste0(celltype,'_lfc'),paste0(celltype,'_FDR'),paste0(celltype,'_flag'))    
    if (length(final)==0){
        final=deseq2
    }else{
        final = merge(final,deseq2)
    }
}
tmp = merge(final,TGF)
dim(tmp)
tmp2 = tmp[order(tmp[,'MN_lfc']),]
data=c()
for (i in seq(1,nrow(tmp2))){
    for (celltype in c('CM','MN','PP','HP','SKM','RPE')){
        coln=paste0(celltype,'_lfc')
        flagn=paste0(celltype,'_FDR')
        if(tmp2[i,flagn]<=0.005){
            labelext = '***'
        }else if(tmp2[i,flagn]<=0.01){
            labelext = '**'
        }else if(tmp2[i,flagn]<=0.05){
            labelext= '*'
        }else{labelext=''}    
        
        data = rbind(data,c(tmp2[i,'gene'],tmp2[i,coln],celltype,labelext))
    }
}
data = data.frame(data,stringsAsFactors =F)
colnames(data)=c('gene','LFC','celltype','DE')
data$LFC=as.numeric(data$LFC)
data$celltype=factor(data$celltype,levels=c('CM','HP','MN','PP','RPE','SKM'))

output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/'
pdf(paste0(output_path,'MN_set_LFC.pdf'),useDingbats=FALSE,width=3.5,height=4)
print(max(data$LFC))
print(min(data$LFC))
ggplot(data, aes(factor(celltype), gene)) + 
  geom_tile(aes(fill=LFC),color='black',height=1,width=1)+
  geom_text(aes(label=DE))+
  scale_fill_gradientn(colours = c("dodgerblue2", "white","indianred1"),limits=c(-2.1,2.1))+
  labs(x = "", y = "", fill = "Log2(Human/Chimp)")+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  theme(legend.position = "top")
dev.off()
