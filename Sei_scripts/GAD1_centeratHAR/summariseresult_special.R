library(stringr)
folder_name='GAD1_centeratHAR'

input_path=paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei/')                          
output_path=paste0('/scratch/users/banwang3/references/outputs/',folder_name,'/sei/')
id='hg38_all'
input_path1=paste0(input_path,id,'/')
df1 = list.files(input_path1,pattern='tsv') 
df1=read.csv(paste0(input_path1,df1),sep='\t')
df1 = data.frame(t(df1)) 
nSNP=ncol(df1)-2 
tmpnames=c('humanfasta','chimpfasta',paste0('humanfasta_withSNP',seq(1,nSNP)))
colnames(df1)=tmpnames

id = 'human2chimp'
input_path1=paste0(input_path,id,'/')
df2 = list.files(input_path1,pattern='tsv')
df2=read.csv(paste0(input_path1,df2),sep='\t')
df2 = data.frame(t(df2))
colnames(df2)='human2chimpfasta'
df = cbind(df1,df2)



ids = c('hg38_all','human2chimp','pantro6_all','chimp2human')
df = c()
for (id in ids){
    print(id)
    input_path1=paste0(input_path,id,'/')
    df1 = list.files(input_path1,pattern='tsv')  
    df1=read.csv(paste0(input_path1,df1),sep='\t')
    df1 = data.frame(t(df1))
    if (id=='hg38_all'){
        nSNP=ncol(df1)-2
        tmpnames=c('humanfasta','chimpfasta',paste0('humanfasta_withSNP',seq(1,nSNP)))
    }else if (id=='human2chimp'){
        tmpnames=c('human2chimpfasta')
    }else if (id =='pantro6_all'){
        nSNP=ncol(df1)
        tmpnames=c(paste0('chimpfasta_withSNP',seq(1,nSNP)))
    }else{
        tmpnames=('chimp2humanfasta')
    }
    colnames(df1)=tmpnames
    if (length(df)==0){
        df = df1
    }else{
        df = cbind(df,df1)
    }
}
                                                                                                                            
df = df[-c(1,2),]
for (i in seq(1,ncol(df))){
    df[,i]=as.numeric(df[,i])
}
tmpnumber=ncol(df)
df$id = row.names(df)
df = df[,c(tmpnumber+1,seq(1,tmpnumber))]
SNPlfc = c()
ref = df[,'humanfasta']
alt = df[,'chimpfasta']
SNPlfc=c(SNPlfc,log2(ref/alt))
for (species in c('human','chimp')){
    idx = grep(paste0(species,'fasta_withSNP'),colnames(df))
    SNPname=colnames(df)[idx[length(idx)]]
    nSNP=as.numeric(str_split_fixed(SNPname,'withSNP',2)[,2])
    if (species=='human'){
       ref = df[,'humanfasta'] 
       for (i in seq(1,nSNP)){
           alt = df[,paste0(species,'fasta_withSNP',i)]
           SNPlfc=cbind(SNPlfc,log2(ref/alt))
        }
       SNPlfc=cbind(SNPlfc,log2(ref/df[,'human2chimpfasta']))
       nSNPhuman=nSNP
    }else{
        alt = df[,'chimpfasta']
        for (i in seq(1,nSNP)){
            ref=df[,paste0(species,'fasta_withSNP',i)]
            SNPlfc=cbind(SNPlfc,log2(ref/alt))
        }
        SNPlfc=cbind(SNPlfc,log2(df[,'chimp2humanfasta']/alt))
        nSNPchimp=nSNP
    }
}
colnames(SNPlfc)=c('hum_chp_fastalfc',paste0('humanbased_humchplfc_withSNP',seq(1,nSNPhuman)),'hum_hum2chp_fastalfc',paste0('chimpbased_humchplfc_withSNP',seq(1,nSNPchimp)),'chp2hum_chp_fastalfc')                
df = cbind(df,SNPlfc)      
write.table(df,file='GAD1_centeratHAR_summary_SEI.txt',sep='\t',quote=F,row.names=F,col.names=T)

