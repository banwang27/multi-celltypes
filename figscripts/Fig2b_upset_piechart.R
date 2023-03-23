library(ggplot2)
library(ComplexUpset)
output_path= '/scratch/users/banwang3/celllines/scripts/RScripts/Rplots10/'
celltypes = c('CM','SKM','HP','PP','MN','RPE')
asegenes<-c()
for (ct in celltypes){
    input_filename=paste0('/scratch/users/banwang3/celllines/analysis/DESeq2_norm/merged/',ct,'_merged_qvalue_DESeq2.txt')
    deseq2 = read.csv(input_filename,header=T,sep='\t')
    asegenes<-c(asegenes,deseq2$gene)
}
asegenes<-unique(asegenes)
df <- data.frame(matrix(ncol = length(celltypes)+1, nrow = length(asegenes)))
colnames(df) <- c('gene',celltypes)
df$gene<-asegenes
df[is.na(df)] <- 0
for (ct in celltypes){ 
    input_filename=paste0('/scratch/users/banwang3/celllines/analysis/DESeq2_norm/merged/',ct,'_merged_qvalue_DESeq2.txt')
    deseq2 = read.csv(input_filename,header=T,sep='\t')  
    colidx = grep(ct, colnames(df))
    df[df$gene %in% deseq2$gene,colidx]=1
}
write.table(df,file='ASE_genes_allcts.txt',sep='\t',quote=F,row.names=F)
pdf(paste0(output_path,'ASE_upset.pdf'), width=10,height=7,useDingbats=FALSE)
upset(
      df,
      celltypes,
      queries=list(
        upset_query(set='CM', fill='#86B875'),
        upset_query(set='HP', fill='#E495A5'),
        upset_query(set='PP', fill='#C7A76C'),
        upset_query(set='MN', fill='#39BEB1'),
        upset_query(set='SKM', fill='#7DB0DD'),
        upset_query(set='RPE', fill='#CD99D8')
        ),
      base_annotations=list(
                        'Intersection Size'=(intersection_size(width=0.5,counts=FALSE)
                                             +scale_y_continuous(expand=expansion(mult=c(0, 0.05)),breaks=seq(0,1500,250))# add some space on the top of the bars
                                             +theme(axis.ticks.y=element_line())
                                             +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour='black'))
                                             )
                        ),
      #stripes=upset_stripes(geom=geom_segment(size=0)),
      matrix=intersection_matrix(geom=geom_point(shape='circle filled',size=2,stroke=0.05)),
      set_sizes=(upset_set_size(geom=geom_bar(width=0.3))
                                + theme(axis.line.x=element_line(colour='black'),axis.ticks.x=element_line(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
                ),
      sort_sets='ascending',
      sort_intersections='descending'
    )
dev.off()

df$sum=rowSums(df[,-1])
data=data.frame(df %>% group_by(sum) %>% summarise(n=n()))
data <- data %>% 
  arrange(desc(sum)) %>%
  mutate(prop = n / sum(data$n) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
data$perc=paste0(signif(data$n/nrow(df)*100,4),'%')  
data$sum=factor(data$sum,levels=seq(1,6)) 
data$lab=paste0(data$sum,'\n(',data$n,')')
pdf(paste0(output_path,'ASE_piechart.pdf'), width=4,height=4,useDingbats=FALSE)   
ggplot(data, aes(x="", y=prop, fill=sum)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() + 
  geom_text(aes(y = ypos, label = lab), color = "black", size=2.8)+ 
  scale_fill_brewer(type="seq",palette=3)+
  theme(legend.position="none")  
dev.off()
