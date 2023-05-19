library(ggplot2)  
library(stringr)
library(dplyr)
library(patchwork)
input_path='/scratch/users/banwang3/celllines/ATAC/human/bam/'
xaxis=seq(170822915,170823954)
celltypes = c('CM','HP','MN','PP','SKM')
for (celltype in celltypes){
  df = read.csv(paste0(input_path,celltype,'_Hume_Chpe_bigwigCompare.bedgraph'),sep='\t',header=F)
  print(max(abs(df$V4)))
  p<-ggplot(data=df,aes(x=V3,y=V4))+
    geom_rect(data=df[1,],aes(xmin=170823007, xmax=170823211, ymin=-Inf, ymax=Inf),fill = alpha("#e4dc2a",0.3))+
    geom_line(color="grey")+
    xlim(170822915, 170823954)+
    ylim(-0.63, 0.63)+
    geom_hline(yintercept=0)+
    xlab('')+
    ylab('')+
    theme_bw()+
    theme(legend.position="none",axis.text.x= element_blank(),axis.ticks.x= element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  pdf(paste0('/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/bw_compare_',celltype,'.pdf'),width=4,height=3, useDingbats=FALSE)
  print(p)
  dev.off()
}