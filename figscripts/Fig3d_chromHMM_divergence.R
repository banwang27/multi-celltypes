library(stringr)
library(ggplot2)
library(dplyr)

celltypes=c('CM','MN','HP','PP','SKM')
df = c()
for (celltype in celltypes){
    peaks = read.csv(paste0('/scratch/users/banwang3/celllines/Enrichment_test/ATAC_fortau/Filtered/',celltype,'_ATAC_Filtered_All_Peaks_fixed.txt'),sep='\t')
    peaks = peaks[,c('Peak.Humreffed','l2fc_Humr','l2fc_Chpr','Human.binom_pval','Max.FDR')]
    peaks$celltype=celltype
    peaks$type=ifelse(grepl('promoter',peaks$Peak.Humreffed),'Promoter','Enhancer')
    df = rbind(df,peaks)
}
ref = read.csv('/scratch/users/banwang3/celllines/ATAC/findSNPs/All_Peaks_Peaklist_in_chromatinstate_group.txt',sep='\t')
#Keep all states, so there could be multiple states in one peak
#Count how many states in each peak
ref= ref[,c('peak_id','Group')]
colnames(ref)=c('Peak.Humreffed','chromatin_state')
ref=ref %>% group_by(Peak.Humreffed,chromatin_state)%>%summarise(n=n())
ref=data.frame(ref)

df1 = merge(df,ref)
df$type=factor(df$type,levels=c('Promoter','Enhancer'))
df1$type = factor(df1$type,levels=c('Promoter','Enhancer'))
df1=df1[df1$chromatin_state!='others',]
write.table(df1,file='All_Peaks_withchromatin_state_info.txt',sep='\t',quote=F,row.names=F,col.names=T)
a = df %>%group_by(celltype,type)%>% summarise(var = var(l2fc_Humr),n=n())
a = a%>% group_by(celltype)%>% mutate(zscore=(var-mean(var))/sd(var))   
b = df1 %>%group_by(celltype,chromatin_state)%>% summarise(var = var(l2fc_Humr),n=n())

b = b%>% filter(n>=2)
b = b%>%group_by(celltype)%>% mutate(zscore=(var-mean(var))/sd(var)) 
df2 = df1[df1$Max.FDR<=0.05,]
b1 = df2%>%group_by(celltype,chromatin_state)%>% summarise(diffCA=n())
b1 = b1%>%filter(diffCA>0)%>%group_by(celltype)%>%mutate(normalized=diffCA/sum(diffCA))

state_annotation=read.csv('/scratch/users/banwang3/celllines/ATAC/findSNPs/state_annotations_processed.csv',sep=',')
group_levels=unique(state_annotation$Group)
group_levels=data.frame(Group=group_levels,paper_name=group_levels)
idx = which(group_levels$Group=='exon')
group_levels$paper_name[idx]='exon and transcription'
idx = which(group_levels$Group=='znf')
group_levels$paper_name[idx]='ZNF genes'
idx = which(group_levels$Group=='promoters')
group_levels$paper_name[idx]='flanking promoter'

result = data.frame(b)
result = merge(result,group_levels,by.x='chromatin_state',by.y='Group')
result1=data.frame(result %>% group_by(paper_name)%>% summarise(median=median(zscore)))
result1=result1[order(result1$median),]
write.table(result1,file='All_Peaks_withchromatin_state_divergent_stat.txt',sep='\t',quote=F,row.names=F,col.names=T)

result$paper_name=factor(result$paper_name,levels=result1$paper_name)

output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/'
pdf(paste0(output_path,'All_Peaks_zscore_chromatin_state_variance.pdf'), useDingbats=FALSE,width=3.5,height=5)
ggplot(result,aes(factor(celltype),paper_name))+
    geom_tile(aes(fill=zscore),color='black')+
    scale_fill_gradientn(colours = c("white","#363680"),guide=guide_colorbar(frame.colour="black",ticks.colour="black", direction="horizontal"))+
    xlab('')+
    ylab('')+
    theme_bw()+
    theme(legend.position = "top")
dev.off()

pdf(paste0(output_path,'All_Peaks_zscore_chromatin_state_variance_flipversion.pdf'), useDingbats=FALSE,width=7,height=3)
ggplot(result,aes(factor(paper_name),celltype))+
    geom_tile(aes(fill=zscore),color='black')+
    scale_fill_gradientn(colours = c("white","#363680"),guide=guide_colorbar(frame.colour="black",ticks.colour="black", direction="vertical"))+
    xlab('')+
    ylab('')+
    scale_y_discrete(limits=rev)+
    theme_bw()+
    theme(legend.position = "right",axis.text.x=element_text(angle=45, hjust=1))
dev.off()
