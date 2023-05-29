library(stringr)
library(ggplot2)
library(patchwork)
input_path='/scratch/users/banwang3/references/outputs/FABP7_centeratpromoter/sei/'
output_path='/scratch/users/banwang3/references/scripts/seifigures/'

df = read.csv(paste0(input_path,'pantro6_centeratpeak_sei_pred_lfc.txt'),sep='\t')
#idx = which(df[,paste0('chimpfasta')]>=0.5| df[,paste0('chimp2humanfasta')]>=0.5|df[,paste0('humanfasta')]>=0.5|df[,paste0('humanfasta_withindel')]>=0.5|df[,paste0('humanfasta_del1G')]>=0.5|df[,paste0('humanfasta_add3bp')]>=0.5)
idx = which(df[,paste0('humanfasta')]>=0.5| df[,paste0('chimpfasta')]>=0.5)
df1 = df[idx,]
print(nrow(df1))
ids = c("hum_chp_fastalfc","chp2hum_chp_fastalfc","indelhum_chp_fastalfc","hum_indelhum_fastalfc","hum_humdel1G_fastalfc","hum_humadd3bp_fastalfc")
#df2 = df1[order(-df1$hum_humdel1G_fastalfc),]
df2 = df1[order(-df1$hum_chp_fastalfc),] 
print(df2[1:10,c('id',ids)])
cor.test(df1$hum_chp_fastalfc,df1$hum_indelhum_fastalfc)$estimate
cor.test(df1$hum_chp_fastalfc,df1$hum_humdel1G_fastalfc)$estimate
cor.test(df1$hum_chp_fastalfc,df1$hum_humadd3bp_fastalfc)$estimate
cor.test(df1$hum_indelhum_fastalfc,df1$hum_humdel1G_fastalfc)$estimate 
cor.test(df1$hum_indelhum_fastalfc,df1$hum_humadd3bp_fastalfc)$estimate  
nSNP=as.numeric(str_split_fixed(colnames(df)[ncol(df)],'withSNP',2)[,2])
#df1 = df
#df = read.csv(paste0(input_path,'fasta_sei_pred_all34SNPsin4096bp_centeratSNP_lfc.txt'),sep='\t')
human=c()
pearsonr=c()
pval=c()
for (i in seq(1,nSNP)){ 
    
    sei = data.frame(lfc = df1[,paste0('chimpbased_humchplfc_withSNP',i)],SNP = paste0('SNP',i),type='chimpbasedSNP')
    
    #corresult=cor.test(df1[,paste0('fastalfc')],df1[,paste0('chimpbased_humchplfc_withSNP',i)])
    #corresult=cor.test(df1[,paste0('chp2hum_chp_fastalfc')],df1[,paste0('chimpbased_humchplfc_withSNP',i)])
    #pearsonr=c(pearsonr,corresult$estimate)
    #pval=c(pval,corresult$p.value)
    #sei2 = data.frame(lfc = df1[,paste0('fastalfc_centeratSNP',i)],SNP = paste0('SNP',i),type='fastalfc') 
    #sei1 = data.frame(lfc = fastasei[,paste0('chimpbased_humchplfc_SNP',i)],SNP = paste0('SNP',i),type='chimpbased',tool='sei')
    #human = rbind(human,deepsea,sei)
    human=rbind(human,sei)
    #chimp=rbind(chimp,deepsea1,sei1)
}
#sei = data.frame(lfc = df1[,'fastalfc'],SNP='Fasta',type='fasta')
#ids = c("hum_chp_fastalfc","chp2hum_chp_fastalfc","indelhum_chp_fastalfc","hum_indelhum_fastalfc")
for (id in ids){
    sei=data.frame(lfc = df1[,id],SNP=id,type='fasta')
    human = rbind(human,sei)
}
#seilfc = log2(fastasei$humanfasta/fastasei$chimpfasta)

#tmp2 = data.frame(lfc = df$fastalfc,SNP='fasta',type='fasta',tool='sei')
#human = rbind(human,tmp2)
#chimp = rbind(chimp,tmp,tmp2)
human$SNP = factor(human$SNP,levels=c(ids,paste0('SNP',seq(1,nSNP))))
#human$SNP = factor(human$SNP,levels=c('Fasta',paste0('SNP',seq(1,nSNP))))
human$type = factor(human$type,levels=c('fasta','chimpbasedSNP'))
#pearsonr=data.frame(pearsonr=pearsonr,pval=pval,SNP=paste0('SNP',seq(1,nSNP)))
#pearsonr$plab<-ifelse(pearsonr$pval<0.001,'***',ifelse(pearsonr$pval<0.01,'**',ifelse(pearsonr$pval<0.05,'*','')))
#pearsonr$SNP=factor(pearsonr$SNP,levels=c(paste0('SNP',seq(1,nSNP))))
#pdf(paste0(output_path,'boxplot_fasta_humanbased_lfc_filtered_all34SNPs_centeratSNP.pdf'),useDingbats=FALSE,width=20,height=3)
#pdf(paste0(output_path,'boxplot_fasta_humanbased_lfc_filtered_all34SNPs_centeratHAR.pdf'),useDingbats=FALSE,width=15,height=4)
p1<-ggplot(human,aes(x = SNP,y=lfc,fill=type))+
     geom_hline(yintercept=0, color = "grey")+
    geom_boxplot(outlier.size = 0.15,show.legend=FALSE)+
     ylab('Log2(Human/Chimp)')+
     xlab('')+
     theme_bw()+
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
#p2<-ggplot(pearsonr,aes(x = SNP,y=pearsonr))+ 
#    geom_hline(yintercept=0, color = "grey")+ 
#    geom_bar(stat = "identity",color="black",fill="white",width=0.5)+                           
#        geom_text(aes(x=SNP,label=plab),vjust="inward")+
#        ylim(-1,1)+
#        ylab('Pearson r')+
#        xlab('')+
#        theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
#p<-p1/p2
p<-p1
 pdf(paste0(output_path,'FABP7_pantro6_centeratpeak_sei_pred_special.pdf'),useDingbats=FALSE,width=15,height=8)
print(p)
dev.off()
 
 data = data.frame(x=human[human$SNP=='Fasta','lfc'],y=human[human$SNP=='SNP21','lfc'])
texts=getCorString(data,'x','y')
texts$x=-Inf
texts$y=Inf
#pdf(paste0(output_path,'scatterplot_fasta_humanbased_lfc_filtered_centeratHAR_SNP21.pdf'),useDingbats=FALSE,width=3,height=3)
pdf(paste0(output_path,'scatterplot_fasta_humanbased_lfc_centeratHAR_SNP21.pdf'),useDingbats=FALSE,width=3,height=3)
ggplot(data,aes(x=x,y=y))+
    geom_smooth(method=lm,se=FALSE, fullrange=TRUE,linewidth=0.5,alpha=0.2,linetype = "longdash")+
    geom_point(size=0.2,alpha=0.3,show.legend = FALSE)+
    geom_text(data=texts,aes(label=text),hjust=0,vjust=1,size=3)+
    xlab('Log2(Human/Chimp) - Fasta')+
    ylab('Log2(Human/Chimp) - SNP21')+
    theme_bw()+
    theme(axis.title = element_text(size = 8),legend.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
 pdf(paste0(output_path,'boxplot_fasta_chimpbased_lfc_filtered.pdf'),useDingbats=FALSE,width=4,height=3)
ggplot(chimp,aes(x = SNP,y=lfc,fill=tool))+
    geom_boxplot(outlier.size = 0.25)+
    ylab('fasta chimpbased pred lfc')+
    xlab('')+ 
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
i=21
idx = which(df[,paste0('humanfasta_centeratSNP',i)] >=0.5| df[,paste0('chimpfasta_centeratSNP',i)]>=0.5)
df1 = df[idx,]
texts=getCorString(df1,'fastalfc_centeratSNP21','humanbased_humchplfc_centeratSNP21') 
texts$fastalfc_centeratSNP21=-Inf
texts$humanbased_humchplfc_centeratSNP21=Inf
pdf(paste0(output_path,'scatterplot_fasta_humanbased_lfc_filtered_SNP21.pdf'),useDingbats=FALSE,width=4,height=3)
ggplot(df1,aes(x=fastalfc_centeratSNP21,y=humanbased_humchplfc_centeratSNP21))+
    geom_smooth(method=lm,se=FALSE, fullrange=TRUE,linewidth=0.5,alpha=0.2,linetype = "longdash")+
    geom_point(size=0.2,alpha=0.3,show.legend = FALSE)+
    geom_text(data=texts,aes(label=text),hjust=0,vjust=1,size=3)+
    xlab('Log2(Human/Chimp) based on fasta')+
    ylab('Log2(Human/Chimp) based on SNP 5')+
    theme_bw()+
    theme(legend.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
idx = grep('SNP21',colnames(df1)) 
idx1 = grep('id',colnames(df1))  
df2 = df1[,c(idx1,idx)]
df2 = df2[order(df2$fastalfc_centeratSNP21),]

data = t(scale(t(as.matrix(df[,-ncol(df)]))))
#data = apply(df, 1, scale)
#sabot=data[order(apply(data, 1, var), decreasing=TRUE),]
sabot = data
#colnames(sabot)=celltypes
rownames(sabot) = str_split_fixed(df$TF,'_',2)[,1]
annotation<-data.frame(cellType=str_split_fixed(colnames(sabot),'_',2)[,1],Allele =str_split_fixed(colnames(sabot),'_',2)[,2])
hc<-as.dist(1 - cor(t(sabot),method='pearson'))
annoCol<-list(cellType=c(CM="#86B875",HP="#E495A5",PP="#C7A76C",MN="#39BEB1",SKM="#7DB0DD"))
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots11/'
pdf(paste0(output_path,'heatmap_pearson_TF_tobias_output_allHyb1.pdf'), useDingbats=FALSE,height=15,width=7)
rownames(annotation)<-colnames(sabot)
out<-pheatmap(as.matrix(sabot),
              show_rownames=T,
              annotation=annotation,
              show_colnames=F,
              annotation_colors = annoCol,
              clustering_distance_rows = hc,
              cutree_rows=10,
              )
dev.off()
