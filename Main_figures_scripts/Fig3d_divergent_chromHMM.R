library(stringr)
library(ggplot2)
library(dplyr)
df = c()

celltypes=c('CM','MN','HP','PP','SKM')
for (celltype in celltypes){
    peaks = read.csv(paste0('/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Filtered/',celltype,'_ATAC_Filtered_All_Peaks_fixed.txt'),sep='\t')
    peaks = peaks[,c('Peak.Humreffed','l2fc_Humr','l2fc_Chpr','Human.binom_pval','Max.FDR')]
    peaks$celltype=celltype
    peaks$type=ifelse(grepl('promoter',peaks$Peak.Humreffed),'Promoter','Enhancer')
    df = rbind(df,peaks)
}
ref = read.csv('/scratch/users/banwang3/celllines/ATAC/findSNPs/Peaks_inFilteredcounts_chromatinstate_info.txt',sep='\t')
#Keep all states, so there could be multiple states in one peak
#ref1 = ref %>%distinct(peak_id,Group,.keep_all = TRUE)
#ref = ref1[,c('peak_id','Group')]


#Count how many states in each peak
ref= ref[,c('Peak.Humreffed','Group')]
colnames(ref)=c('Peak.Humreffed','chromatin_state')
ref=ref %>% group_by(Peak.Humreffed,chromatin_state)%>%summarise(n=n())
ref=data.frame(ref)
row.has.na <- apply(ref, 1, function(x){any(is.na(x))})
ref=ref[!row.has.na,]

df1 = merge(df,ref)
df1=df1[df1$chromatin_state!='others',]
write.table(df1,file='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Filtered/All_Peaks_withchromatin_state_info.txt',sep='\t',quote=F,row.names=F,col.names=T)

state_annotation=read.csv('/scratch/users/banwang3/celllines/ATAC/findSNPs/state_annotations_processed.csv',sep=',')
group_levels=unique(state_annotation$Group)
group_levels=data.frame(Group=group_levels,paper_name=group_levels)
df1 = merge(df1,group_levels,by.x='chromatin_state',by.y='Group') 
idx = which(group_levels$Group=='exon')
group_levels$paper_name[idx]='exon and transcription'
idx = which(group_levels$Group=='znf')
group_levels$paper_name[idx]='ZNF genes'
idx = which(group_levels$Group=='promoters')
group_levels$paper_name[idx]='flanking promoter'

df1 = read.csv('/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Filtered/All_Peaks_withchromatin_state_info.txt',sep='\t')
b = df1 %>%group_by(celltype,chromatin_state)%>% summarise(var = var(l2fc_Humr),absmedian=median(abs(l2fc_Humr)),absmean=mean(abs(l2fc_Humr)),n=n())
b = b%>%group_by(celltype)%>% mutate(zscore=(absmedian-mean(absmedian))/sd(absmedian)) 
result = data.frame(b)
result = merge(result,group_levels,by.x='chromatin_state',by.y='Group')
idx = which(result$paper_name=='quescient')
result$paper_name[idx]='quiescent'
idx = which(result$paper_name=='transcribed and enhancer')
result$paper_name[idx]='transcribed & enhancer'
idx = which(result$paper_name=='exon and transcription')
result$paper_name[idx]='exon & transcription'
result1=data.frame(result %>% group_by(paper_name)%>% summarise(median=median(zscore)))
result1=result1[order(result1$median),]
write.table(result1,file='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Filtered/All_Peaks_withchromatin_state_divergent_stat.txt',sep='\t',quote=F,row.names=F,col.names=T)


result$paper_name=factor(result$paper_name,levels=result1$paper_name)
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/'
pdf(paste0(output_path,'All_Peaks_chromatin_state_median_absLFC.pdf'), useDingbats=FALSE,width=3.5,height=5)
ggplot(result,aes(factor(celltype),paper_name))+
    geom_tile(aes(fill=zscore),color='black')+
    scale_fill_gradientn(colours = c("white","#363680"),guide=guide_colorbar(frame.colour="black",ticks.colour="black", direction="horizontal"))+
    labs(x = "", y = "", fill = "median|LFC|\nz-score") +
    theme_bw()+
    theme(legend.position = "top")
dev.off()


