folder_name='FABP7_centeratpromoter'
input_path = paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei_pantro6/sei/')
input_path1=paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei_pantro6/sei1/')
input_path2 = paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei/')
input_path_del1G = paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei_del1G/')
input_path_add3bp= paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei_add3bp/')
                          
output_path=paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei/')
align2='pantro6'
df1 = list.files(input_path1,pattern='tsv')
df1=read.csv(paste0(input_path1,'/',df1),sep='\t')
df1 = data.frame(t(df1))
df = list.files(input_path,pattern='tsv')
df= read.csv(paste0(input_path,'/',df),sep='\t')
df = data.frame(t(df)) 
nSNP=ncol(df)-2
df2 = list.files(input_path2,pattern='tsv')
df2=read.csv(paste0(input_path2,'/',df2),sep='\t')
df2= data.frame(t(df2)) 
df_del1G = list.files(input_path_del1G,pattern='tsv') 
df_del1G=read.csv(paste0(input_path_del1G,'/',df_del1G),sep='\t')
df_del1G=data.frame(t(df_del1G))
df_add3bp = list.files(input_path_add3bp,pattern='tsv')  
df_add3bp=read.csv(paste0(input_path_add3bp,'/',df_add3bp),sep='\t') 
df_add3bp= data.frame(t(df_add3bp))
df = cbind(df,df1,df2,df_del1G,df_add3bp)
colnames(df)=c(paste0('chimpfasta'),paste0('humanfasta'),paste0('chimpfasta_withSNP',seq(1:nSNP)),paste0('chimp2humanfasta'),paste0('humanfasta_withindel'),paste0('humanfasta_del1G'),paste0('humanfasta_add3bp'))
                                                                                                                            
df = df[-c(1,2),]
for (i in seq(1,ncol(df))){
    df[,i]=as.numeric(df[,i])
}
df$id = row.names(df)
SNPlfc = c()
ref = df[,'humanfasta']
alt = df[,'chimpfasta']
c2h_ref = df[,'chimp2humanfasta']
withindel_ref=df[,'humanfasta_withindel']
fastalfc=log2(ref/alt)
fastalfc2=log2(c2h_ref/alt)
fastalfc3=log2(withindel_ref/alt)
fastalfc4=log2(ref/withindel_ref)
fastalfc5=log2(ref/df[,'humanfasta_del1G'])
fastalfc6=log2(ref/df[,'humanfasta_add3bp'])
for (i in seq(1,nSNP)){
    altSNP = df[,2+i]
    SNPlfc = cbind(SNPlfc,log2(altSNP/alt))
}
tmpnames=colnames(df)
df = cbind(df,fastalfc,fastalfc2,fastalfc3,fastalfc4,fastalfc5,fastalfc6,SNPlfc)   
colnames(df)=c(tmpnames,paste0('hum_chp_fastalfc'),paste0('chp2hum_chp_fastalfc'),paste0('indelhum_chp_fastalfc'),paste0('hum_indelhum_fastalfc'),paste0('hum_humdel1G_fastalfc'),paste0('hum_humadd3bp_fastalfc'),paste0('chimpbased_humchplfc_withSNP',seq(1:nSNP)))
 write.table(df,file=paste0(output_path,align2,'_centeratpeak_sei_pred_lfc.txt'),sep='\t',quote=F,row.names=F,col.names=T)
