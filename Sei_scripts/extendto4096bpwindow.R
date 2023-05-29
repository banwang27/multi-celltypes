library(argparse)
parser <- ArgumentParser()
parser$add_argument("--generator", default="FABP7")
args <- parser$parse_args()
foldername=args$generator
input_path=paste0('/scratch/users/banwang3/references/outputs/',foldername,'/')
align2='hg38'
for (align2 in c('hg38','pantro6')){
    df = read.csv(paste0(input_path,align2,'_peak_coord.bed'),sep='\t',header=F)
    peakid=df$V4[1]
    if (grepl('promoter',peakid)){
        type='promoter'
    }else{
        type='enhancer'
    }
    peakchr=df$V1[1]
    peakstart = df$V2[1]
    peakend = df$V3[1]
    newstart=floor((peakend-peakstart)/2)+peakstart-2048
    newend=floor((peakend-peakstart)/2)+peakstart+2048
    newpeak = data.frame(chr=peakchr,start=newstart,end=newend,id=paste0(foldername,'_',type))
    write.table(newpeak,file=paste0(input_path,align2,'_centeratpeak_4096bp.bed'),sep='\t',quote=F,row.names=F,col.names=F)
}
