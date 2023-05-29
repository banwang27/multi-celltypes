folder_name='FABP7'
input_path = paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei/')
input_path1=paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei1/')
output_path=paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei/')
align2='hg38'
df1 = list.files(input_path1,pattern='tsv')
df1=read.csv(paste0(input_path1,'/',df1),sep='\t')
df1 = data.frame(t(df1))
df = list.files(input_path,pattern='tsv')
df= read.csv(paste0(input_path,'/',df),sep='\t')
df = data.frame(t(df)) 
nSNP=ncol(df)-2
df = cbind(df,df1)
colnames(df)=c(paste0('humanfasta'),paste0('chimpfasta'),paste0('humanfasta_withSNP',seq(1:nSNP)),paste0('human2chimpfasta'))
df = df[-c(1,2),]
for (i in seq(1,ncol(df))){
    df[,i]=as.numeric(df[,i])
}
df$id = row.names(df)
SNPlfc = c()
ref = df[,1]
alt = df[,2]
h2c_alt = df[,'human2chimpfasta']
fastalfc=log2(ref/alt)
fastalfc2=log2(ref/h2c_alt)
for (i in seq(1,nSNP)){
    refSNP = df[,2+i]
    SNPlfc = cbind(SNPlfc,log2(ref/refSNP))
}
tmpnames=colnames(df)
df = cbind(df,fastalfc,fastalfc2,SNPlfc)   
colnames(df)=c(tmpnames,paste0('fastalfc'),paste0('hum_hum2chp_fastalfc'),paste0('humanbased_humchplfc_withSNP',seq(1:nSNP)))
 write.table(df,file=paste0(output_path,align2,'_centeratpeak_sei_pred_lfc.txt'),sep='\t',quote=F,row.names=F,col.names=T)
