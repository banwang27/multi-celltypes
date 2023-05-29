folder_name='FABP7_centeratpromoter'
input_path = paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei_pantro6/sei/')
input_path1=paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei_pantro6/sei1/')
output_path=paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei_pantro6/sei/')
align2='pantro6'
df1 = list.files(input_path1,pattern='tsv')
df1=read.csv(paste0(input_path1,'/',df1),sep='\t')
df1 = data.frame(t(df1))
df = list.files(input_path,pattern='tsv')
df= read.csv(paste0(input_path,'/',df),sep='\t')
df = data.frame(t(df)) 
nSNP=ncol(df)-2
df = cbind(df,df1)
colnames(df)=c(paste0('chimpfasta'),paste0('humanfasta'),paste0('chimpfasta_withSNP',seq(1:nSNP)),paste0('chimp2humanfasta'))
df = df[-c(1,2),]
for (i in seq(1,ncol(df))){
    df[,i]=as.numeric(df[,i])
}
df$id = row.names(df)
SNPlfc = c()
ref = df[,'humanfasta']
alt = df[,'chimpfasta']
c2h_ref = df[,'chimp2humanfasta']
fastalfc=log2(ref/alt)
fastalfc2=log2(c2h_ref/alt)
for (i in seq(1,nSNP)){
    altSNP = df[,2+i]
    SNPlfc = cbind(SNPlfc,log2(altSNP/alt))
}
tmpnames=colnames(df)
df = cbind(df,fastalfc,fastalfc2,SNPlfc)   
colnames(df)=c(tmpnames,paste0('fastalfc'),paste0('chp2hum_chp_fastalfc'),paste0('chimpbased_humchplfc_withSNP',seq(1:nSNP)))
 write.table(df,file=paste0(output_path,align2,'_centeratpeak_sei_pred_lfc.txt'),sep='\t',quote=F,row.names=F,col.names=T)
