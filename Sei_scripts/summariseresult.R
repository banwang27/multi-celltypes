library(argparse)
parser <- ArgumentParser()
parser$add_argument("--generator", default="GLI2_1")
args <- parser$parse_args()
folder_name=args$generator
input_path = paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei/')
output_path=paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei/')
align2='hg38'
filename=paste0(input_path,align2,'_all')
df = list.files(filename,pattern='tsv')
df= read.csv(paste0(filename,'/',df),sep='\t')
df = data.frame(t(df))
nSNP=ncol(df)/3
colnames(df)=c(paste0('humanfasta_centeratSNP',seq(1:nSNP)),paste0('chimpfasta_centeratSNP',seq(1:nSNP)),paste0('humanfasta_withSNP',seq(1:nSNP)))
df = df[-c(1,2),]
for (i in seq(1,ncol(df))){
    df[,i]=as.numeric(df[,i])
}
df$id = row.names(df)
fastalfc = c()
SNPlfc = c()
for (i in seq(1,nSNP)){
    ref = df[,i]
    alt = df[,i+nSNP]
    refSNP = df[,i+nSNP+nSNP]
    fastalfc = cbind(fastalfc,log2(ref/alt))
    SNPlfc = cbind(SNPlfc,log2(ref/refSNP))
}
tmpnames=colnames(df)
df = cbind(df,fastalfc,SNPlfc)   
colnames(df)=c(tmpnames,paste0('fastalfc_centeratSNP',seq(1:nSNP)),paste0('humanbased_humchplfc_centeratSNP',seq(1:nSNP)))
 write.table(df,file=paste0(output_path,align2,'_centeratSNP_sei_pred_lfc.txt'),sep='\t',quote=F,row.names=F,col.names=T)
