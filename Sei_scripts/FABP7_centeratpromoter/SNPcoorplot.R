library(ggplot2)
library(scales)   
output_path='/scratch/users/banwang3/references/scripts/seifigures/'
SNP=read.csv('SNPlist_hg38_4096bp.bed',sep='\t',header=F)
windowstart = SNP$V6[1]
windowend= SNP$V7[1]
SNP=SNP[,1:3]
colnames(SNP)=c('chr','start','end')
SNP$label=paste0('SNP',seq(1,nrow(SNP)))
SNP$order=seq(1,nrow(SNP))
SNP$y=ifelse(SNP$order%%2==0,0.4,0.1)
highlight=read.csv('hg38_peak_coord.bed',sep='\t',header=F)
highlightstart=highlight$V2[1]
highlightend=highlight$V3[1]

p<-ggplot(data=SNP)+
    geom_rect(data=SNP[1,],aes(xmin=highlightstart,xmax=highlightend, ymin=-Inf, ymax=Inf),fill = alpha("#e4dc2a",0.3))+
    geom_segment( aes(x=end, xend=end, y=0, yend=y),color="grey")+
    geom_point(data=SNP,aes(x=end,y=y),color="red",size=2)+
    geom_hline(yintercept=0, color = "black")+
    xlab('')+
    ylim(0,0.5)+
    scale_x_continuous(limits=c(windowstart, windowend), expand = c(0, 0),labels=comma,breaks=seq((floor(windowstart/500)+1)*500,floor(windowend/500)*500,500)) +
    theme_bw()+ 
    theme(panel.border = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())    
pdf(paste0(output_path,'FABP7_SNP_coor.pdf'),width=9.5,height=0.8, useDingbats=FALSE)
print(p)
dev.off()

