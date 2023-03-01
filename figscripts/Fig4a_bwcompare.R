library(ggplot2)  
library(stringr)
library(dplyr)
library(patchwork)
input_path='/scratch/users/banwang3/celllines/ATAC/human/bam/'
HARstart = 170823007
HARend = 170823211
peakstart = 170822915
peakend = 170823954
celltypes = c('CM','HP','MN','PP','SKM')
for (celltype in celltypes){
    df = read.csv(paste0(input_path,celltype,'_Hume_Chpe_bigwigCompare.bedgraph'),sep='\t',header=F)
    print(max(abs(df$V4)))
    p<-ggplot(data=df,aes(x=V3,y=V4))+
        geom_rect(data=df[1,],aes(xmin=HARstart, xmax=HARend, ymin=-Inf, ymax=Inf),fill = alpha("#e4dc2a",0.3))+
        geom_line(color="grey")+
        xlim(peakstart, peakend)+
        ylim(-0.63, 0.63)+
        geom_hline(yintercept=0)+
        xlab('')+
        #ylab('Log(Human/Chimp)')+
        ylab('')+
        theme_bw()+
        theme(legend.position="none",axis.text.x= element_blank(),axis.ticks.x= element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    pdf(paste0('/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/bw_compare_',celltype,'.pdf'),width=4,height=3, useDingbats=FALSE)
    print(p)
    dev.off()       
}
