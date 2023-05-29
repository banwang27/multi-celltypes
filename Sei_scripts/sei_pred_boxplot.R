library(stringr)
library(ggplot2)
library(patchwork)
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--generator", default="FABP7")
parser$add_argument("--cutoff",type="integer", default=0.5)
args <- parser$parse_args()
foldername=args$generator
cutoff=args$cutoff
input_path=paste0('/scratch/users/banwang3/references/outputs/',foldername,'/sei/')
output_path='/scratch/users/banwang3/references/scripts/seifigures/'
align2='hg38'
df = read.csv(paste0(input_path,align2,'_centeratSNP_sei_pred_lfc.txt'),sep='\t')
nSNP=as.numeric(str_split_fixed(colnames(df)[ncol(df)],'SNP',2)[,2]) 
human=c()
pearsonr=c()
pval=c()
count=c()
for (i in seq(1,nSNP)){ 
    idx = which(df[,paste0('humanfasta_centeratSNP',i)] >=cutoff| df[,paste0('chimpfasta_centeratSNP',i)]>=cutoff)
    df1 = df[idx,] 
    count=c(count,nrow(df1))
    if (nrow(df1)<=2){
            fastalfc=data.frame(lfc=0,SNP=paste0('SNP',i),type='fasta',count=nrow(df1))
            sei = data.frame(lfc=0,SNP = paste0('SNP',i),type='humanbasedSNP',count=nrow(df1))
          pearsonr=c(pearsonr,0)
          pval=c(pval,1)
     }else{
         fastalfc=data.frame(lfc=df1[,paste0('fastalfc_centeratSNP',i)],SNP=paste0('SNP',i),type='fasta',count=nrow(df1))
         sei = data.frame(lfc = df1[,paste0('humanbased_humchplfc_centeratSNP',i)],SNP = paste0('SNP',i),type='humanbasedSNP',count=nrow(df1))
          corresult=cor.test(df1[,paste0('fastalfc_centeratSNP',i)],df1[,paste0('humanbased_humchplfc_centeratSNP',i)])
          pearsonr=c(pearsonr,corresult$estimate)
         pval=c(pval,corresult$p.value)
     }
    human=rbind(human,fastalfc,sei)
}
human$SNP = factor(human$SNP,levels=c(paste0('SNP',seq(1,nSNP))))
human$type = factor(human$type,levels=c('fasta','humanbasedSNP'))
pearsonr=data.frame(pearsonr=pearsonr,pval=pval,SNP=paste0('SNP',seq(1,nSNP)))
pearsonr$plab<-ifelse(pearsonr$pval<0.001,'***',ifelse(pearsonr$pval<0.01,'**',ifelse(pearsonr$pval<0.05,'*','')))
pearsonr$SNP=factor(pearsonr$SNP,levels=c(paste0('SNP',seq(1,nSNP))))
p1<-ggplot(human,aes(x = SNP,y=lfc,fill=type))+
     geom_hline(yintercept=0, color = "grey")+
    geom_boxplot(outlier.size = 0.15,show.legend=T)+
     ylab('Log2(Human/Chimp)')+
     xlab('')+
     theme_bw()+
     theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p2<-ggplot(pearsonr,aes(x = SNP,y=pearsonr))+                                             
      geom_hline(yintercept=0, color = "grey")+                                         
      geom_bar(stat = "identity",color="black",fill="white",width=0.5)+                           
      geom_text(aes(x=SNP,label=plab),vjust="inward")+
      ylim(-1,1)+
      ylab('Pearson r')+
      xlab('')+
      theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p<-p1/p2 
pdf(paste0(output_path,foldername,'_',align2,'_sei_centeratSNP.pdf'),useDingbats=FALSE,width=15,height=8)
print(p) 
dev.off()

