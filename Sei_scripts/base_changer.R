library(stringr)
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--generator", default="FABP7")
args <- parser$parse_args()
foldername=args$generator
input_path=paste0('/scratch/users/banwang3/references/outputs/',foldername,'/')
align2s=c('hg38','pantro6')
for (align2 in align2s){
filename=paste0('SNPlist_',align2,'_4096bp.bed')
SNPs=read.csv(paste0(input_path,filename),sep='\t',header=F)
SNPs=SNPs[,1:4]
colnames(SNPs)=c('chr','start','end','info')
SNPs$ref = str_split_fixed(SNPs$info,"\\|",2)[,1]
SNPs$alt = str_split_fixed(SNPs$info,"\\|",2)[,2]
print(SNPs[1,])
fasta = readLines(paste0(input_path,align2,'_centeratSNP_4096bp.fasta'))
count=0
towrite=c()
for (i in seq(1,nrow(SNPs))){
    fastaseq = fasta[i*2]
    pos = substring(fastaseq,2048,2048)
    
    newstring = fastaseq
    if (align2=='hg38'){
        baseallele=SNPs$ref[i]
        changeallele=SNPs$alt[i]
    }else{
        baseallele=SNPs$alt[i]
        changeallele=SNPs$ref[i]
    }
    if (toupper(pos)==baseallele){
        substring(newstring,2048,2048)<-changeallele
        count=count+1
    }
    towrite=c(towrite,paste0('>SNPfasta_centeratSNP',i,'_align2',align2),newstring)
}
if (count==nrow(SNPs)){
    print('all SNPs written')
}
filename=paste0(input_path,'SNPfasta_',align2,'_centeratSNP_4096bp.fasta')
fileConn<-file(filename)
writeLines(towrite, fileConn)
close(fileConn)
}
