library(stringr)
#library(dplyr)
foldername='GAD1_centeratHAR'
input_path=paste0('/scratch/users/banwang3/references/outputs/',foldername,'/')
align2='hg38'
#align2s=c('hg38','pantro6')
filename=paste0('SNPlist_',align2,'_4096bp.bed')
SNPs=read.csv(paste0(input_path,filename),sep='\t',header=F)
SNPs=SNPs[,1:4]
colnames(SNPs)=c('chr','start','end','info')
SNPs$ref = str_split_fixed(SNPs$info,"\\|",2)[,1]
SNPs$alt = str_split_fixed(SNPs$info,"\\|",2)[,2]
print(SNPs[1,])
fasta = readLines(paste0(input_path,align2,'_centeratHAR_4096bp.fasta'))
fasta1 = readLines(paste0(input_path,'pantro6','_centeratHAR_4096bp.fasta'))
fastainfofull=fasta[1]
fastaseq = fasta[2]
fastainfo=str_split_fixed(fastainfofull,':',2)[,2]
fastastart = as.numeric(str_split_fixed(fastainfo,'-',2)[,1])
fastaend= as.numeric(str_split_fixed(fastainfo,'-',2)[,2])
fastainfo1=fasta1[1]
fastaseq1=fasta1[2]
df = data.frame(idx=seq(1,4096),human=strsplit(fastaseq, "")[[1]],chimp=strsplit(fastaseq1, "")[[1]])
df$human=toupper(df$human)
df$chimp=toupper(df$chimp)
2971
df = df[2971:nrow(df),]
df = data.frame(idx=seq(1:(nrow(df)-1)),human=df$human[2:nrow(df)],chimp=df$chimp[1:(nrow(df)-1)])
574
df = df[574:nrow(df),]
df = data.frame(idx=seq(1:(nrow(df)-3)),human=df$human[4:nrow(df)],chimp=df$chimp[1:(nrow(df)-3)])
idx = which(df$human!=df$chimp)
#fasta2=readLines(paste0(input_path,'human2chimp','_centeratHAR_4096bp.fasta'))
#fastainfo2=fasta2[1]
#fastaseq2=fasta2[2]
#towrite=c(fastainfofull,fastaseq,fastainfo1,fastaseq1,fastainfo2,fastaseq2)
towrite=c(fastainfofull,fastaseq,fastainfo1,fastaseq1)
count=0
towrite=c()
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
    towrite=c(towrite,paste0('>SNPfasta_centeratHAR',i,'_align2',align2),newstring)
    if (count==nrow(SNPs)){
        print('all SNPs written')
    }
}
filename=paste0(input_path,'SNPfasta_',align2,'_centeratHAR_4096bp.fasta')
fileConn<-file(filename)
writeLines(towrite, fileConn)
close(fileConn)

