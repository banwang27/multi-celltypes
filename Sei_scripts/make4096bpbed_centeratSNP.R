library(argparse)
parser <- ArgumentParser()
parser$add_argument("--generator", default="FABP7")
args <- parser$parse_args()
foldername=args$generator
input_path=paste0('/scratch/users/banwang3/references/outputs/',foldername,'/')
for (align2 in c('hg38','pantro6')){
    filename = paste0(input_path,'SNPlist_',align2,'_4096bp.bed')
    SNP=read.csv(filename,sep='\t',header=F)
    print(SNP[1,])
    SNP=SNP[,c(1,2,3,4)]
    colnames(SNP)=c('chr','start','end','info')
    newbed = cbind(SNP$chr,SNP$end-2048,SNP$end+2048,paste0(align2,'_centeratSNP',seq(1,nrow(SNP))))
    SNP$count=paste0('SNP',seq(1,nrow(SNP))) 
    peakinfo=read.csv(paste0(input_path,align2,'_peak_coord.bed'),sep='\t',header=F)
    SNPsub=SNP[SNP$start>=peakinfo$V2[1] & SNP$end<=peakinfo$V3[1],]
    write.table(SNPsub,file=paste0(input_path,align2,'_SNPinpeak.bed'),quote=F,sep='\t',row.names=F,col.names=F)
    write.table(newbed,file=paste0(input_path,align2,'_centeratSNP_4096bp.bed'),quote=F,sep='\t',row.names=F,col.names=F)
}


# For lift over version
filename=paste0(input_path,'SNPlist_hg38_4096bp_lift2pantro6.bed')
SNP=read.csv(filename,sep='\t',header=F)
print(SNP[1,])
print(nrow(SNP))
SNP=SNP[,c(1,2,3,4)]
colnames(SNP)=c('chr','start','end','info')
newbed = cbind(SNP$chr,SNP$end-2048,SNP$end+2048,paste0('hg38_lift2pantro6_centeratSNP',seq(1,nrow(SNP))))
write.table(newbed,file=paste0(input_path,'hg38_centeratSNP_4096bp_lift2pantro6.bed'),quote=F,sep='\t',row.names=F,col.names=F)

filename=paste0(input_path,'SNPlist_pantro6_4096bp_lift2hg38.bed')
SNP=read.csv(filename,sep='\t',header=F)
print(SNP[1,])    
print(nrow(SNP))
SNP=SNP[,c(1,2,3,4)] 
colnames(SNP)=c('chr','start','end','info')
newbed = cbind(SNP$chr,SNP$end-2048,SNP$end+2048,paste0('pantro6_lift2hg38_centeratSNP',seq(1,nrow(SNP))))
write.table(newbed,file=paste0(input_path,'pantro6_centeratSNP_4096bp_lift2hg38.bed'),quote=F,sep='\t',row.names=F,col.names=F)
