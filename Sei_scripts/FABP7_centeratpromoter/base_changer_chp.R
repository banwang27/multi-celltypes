library(stringr)
#library(dplyr)
foldername='FABP7_centeratpromoter'
input_path=paste0('/scratch/users/banwang3/references/outputs/',foldername,'/')
align2='pantro6'
#align2s=c('hg38','pantro6')
filename=paste0('SNPlist_',align2,'_4096bp.bed')
SNPs=read.csv(paste0(input_path,filename),sep='\t',header=F)
SNPs=SNPs[,1:4]
colnames(SNPs)=c('chr','start','end','info')
SNPs$ref = str_split_fixed(SNPs$info,"\\|",2)[,1]
SNPs$alt = str_split_fixed(SNPs$info,"\\|",2)[,2]
print(SNPs[1,])
fasta = readLines(paste0(input_path,align2,'_centeratpeak_4096bp.fasta'))
fasta1 = readLines(paste0(input_path,'hg38','_centeratpeak_4096bp.fasta'))
fastainfofull=fasta[1]
fastaseq = fasta[2]
a = unlist(strsplit(fastaseq, split = ""))
fastainfo=str_split_fixed(fastainfofull,':',2)[,2]
fastastart = as.numeric(str_split_fixed(fastainfo,'-',2)[,1])
fastaend= as.numeric(str_split_fixed(fastainfo,'-',2)[,2])
fastainfo1=fasta1[1]
fastaseq1=fasta1[2]
fasta2=readLines(paste0(input_path,'chimp2human','_centeratpeak_4096bp.fasta'))
fastainfo2=fasta2[1]
fastaseq2=fasta2[2]
b = unlist(strsplit(fastaseq2, split = ""))
#towrite=c(fastainfofull,fastaseq,fastainfo1,fastaseq1,fastainfo2,fastaseq2)
towrite=c(fastainfofull,fastaseq,fastainfo1,fastaseq1)
count=0
#towrite=c()
for (i in seq(1,nrow(SNPs))){
    pos = substring(fastaseq,SNPs$end[i]-fastastart,SNPs$end[i]-fastastart)
    
    newstring = fastaseq
    if (align2=='hg38'){
        baseallele=SNPs$ref[i]
        changeallele=SNPs$alt[i]
    }else{
        baseallele=SNPs$alt[i]
        changeallele=SNPs$ref[i]
    }
    if (toupper(pos)==baseallele){
        substring(newstring,SNPs$end[i]-fastastart,SNPs$end[i]-fastastart)<-changeallele
        count=count+1
    }
    towrite=c(towrite,paste0('>SNPfasta_centeratpeak',i,'_align2',align2),newstring)
    if (count==nrow(SNPs)){
        print('all SNPs written')
    }
}
filename=paste0(input_path,'SNPfasta_',align2,'_centeratpeak_4096bp_2.fasta')
fileConn<-file(filename)
writeLines(towrite, fileConn)
close(fileConn)

