#library(UpSetR)
library(ggplot2)
library(ComplexUpset)
output_path= '/scratch/users/banwang3/celllines/scripts/RScripts/suppplots/'
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/'
celltypes = c('CM','SKM','HP','PP','MN')
peaks=c()
for (celltype in celltypes){
    df = read.csv(paste0(input_path,'Downsampled/',celltype,'_ATAC_Filtered_Down_fixed.txt'),sep='\t',header=T)
    df = df[!grepl('chr20',df$Genomic.location),]
    peaks = c(peaks,df$Peak.Humreffed)
}
peaks=unique(peaks)
df1 <- data.frame(matrix(ncol = length(celltypes)+1, nrow = length(peaks)))
colnames(df1)<-c('peak_id',celltypes)
df1$peak_id<-peaks
for (celltype in celltypes){
    df1[,celltype]=ifelse(grepl(celltype,df1$peak_id),1,0)
}
pdf(paste0(output_path,'Allpeaks_downsampled_upset.pdf'), width=10,height=7,useDingbats=FALSE)
upset(
      df1,
      celltypes,
      queries=list(
        upset_query(set='CM', fill='#86B875'),
        upset_query(set='HP', fill='#E495A5'),
        upset_query(set='PP', fill='#C7A76C'),
        upset_query(set='MN', fill='#39BEB1'),
        upset_query(set='SKM', fill='#7DB0DD')
        ),
      base_annotations=list(
                        'Intersection Size'=(intersection_size(width=0.5,counts=FALSE)
                                             #+scale_y_continuous(expand=expansion(mult=c(0, 0.05)),breaks=seq(0,1500,250))# add some space on the top of the bars
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


