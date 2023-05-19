library(stringr)
library(ggplot2)
library(dplyr)
celltypes=c('CM','MN','HP','PP','SKM')
reffed='humreffed'
align2='human'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/'
backup_path='/scratch/users/banwang3/celllines/data_backup'
df = read.csv('/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Filtered/All_Peaks_withchromatin_state_info.txt',sep='\t')
bnomial_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/normalized/Normalized_binom/'
    
final=c()
for (celltype in celltypes){
    dfopen = df[df$celltype==celltype,]
    dff = read.csv(paste0(bnomial_path,celltype,'_Humreffed_NotDown_Filtered_Normed_Binom.txt'),sep='\t')
    dff = dff[,c('Peak.Humreffed','Binomial.p.value')]
    dfopen = merge(dfopen,dff)    
    dfopen$flag = ifelse(dfopen$Binomial.p.value<0.05,'Allele-Specific Peaks','All Peaks')
    dfopen_prom = dfopen[dfopen$type=='Promoter',]
    dfopen_enha = dfopen[dfopen$type=='Enhancer',]
    dfopen_prom$gene = str_split_fixed(dfopen_prom$Peak.Humreffed,'promoter_',2)[,2]
    # there could be one promoter peak targeting multiple genes
    while(length(which(grepl('_',dfopen_prom$gene)))>0){
        idx = which(grepl('_',dfopen_prom$gene))
        uniquegenedf = dfopen_prom[-idx,]
        dfopen_prom = dfopen_prom[idx,]
        gene1 = str_split_fixed(dfopen_prom$gene,'_',2)[,1] 
        gene2 = str_split_fixed(dfopen_prom$gene,'_',2)[,2] 
        tmp1 = dfopen_prom
        tmp1$gene = gene1
        tmp2 = dfopen_prom
        tmp2$gene = gene2
        dfopen_prom = rbind(uniquegenedf,tmp1,tmp2)
    } 
    #For enhancer peaks, take only the closet gene as target gene 
    dfopen_enha$gene = str_split_fixed(dfopen_enha$Peak.Humreffed,'enhancer_',2)[,2]
    dfopen_enha$gene = str_split_fixed(dfopen_enha$gene,'_',2)[,1]
    # get gene ASE value computed by DESeq2
    deseq2 = read.csv(paste0(backup_path,'/scratch/users/banwang3/celllines/analysis/DESeq2_norm/',align2,'/',celltype,'_aligned2human_qvalue_DESeq2.txt'),sep='\t')
    dfopen_prom = merge(dfopen_prom,deseq2,all.x=TRUE)
    dfopen_enha = merge(dfopen_enha,deseq2,all.x=TRUE)
    row.has.na <- apply(dfopen_prom, 1, function(x){any(is.na(x))})
    dfopen_prom = dfopen_prom[!row.has.na,]
    row.has.na <- apply(dfopen_enha, 1, function(x){any(is.na(x))})    
    dfopen_enha = dfopen_enha[!row.has.na,]
    df1 = rbind(dfopen_prom,dfopen_enha)
    ASCA=df1[df1$flag=='Allele-Specific Peaks',]
    df2 = data.frame(ASCA%>% group_by(paper_name)%>%filter(n() >=15)%>%summarise(pearsonr=cor.test(log2FoldChange,l2fc_Humr)$estimate,pval=cor.test(log2FoldChange,l2fc_Humr)$p.val,n=n())) 
    df2$flag='Allele-Specific Peaks'
    df3 = data.frame(df1%>% group_by(paper_name)%>%filter(n() >=15)%>%summarise(pearsonr=cor.test(log2FoldChange,l2fc_Humr)$estimate,pval=cor.test(log2FoldChange,l2fc_Humr)$p.val,n=n()))
    df3$flag='All Peaks'
    df4 = rbind(df2,df3)
    df4$celltype=celltype
    final=rbind(final,df4)
}
statenames=read.csv('/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/ATAC_counttable/Filtered/All_Peaks_withchromatin_state_divergent_stat.txt',sep='\t')
final$labelest=ifelse(final$pval<0.005,'***',ifelse(final$pval<0.01,'**',ifelse(final$pval<0.05,'*','')))
final$paper_name=factor(final$paper_name,levels=c(statenames$paper_name))
print(max(abs(final$pearsonr)))


finalASCA=final[final$flag=='Allele-Specific Peaks',]
print(min(finalASCA$n))
pdf(paste0(output_path,'ASCA_bnomialp_Peaks_chromatin_state_ASEcor_orderbydivergent.pdf'), width=3.5,height=5,useDingbats=FALSE)
ggplot(finalASCA,aes(factor(celltype),paper_name))+
    geom_tile(aes(fill=pearsonr),colour = "black")+
    geom_text(aes(label=labelest),size=2)+
    scale_fill_gradientn(colours = c("dodgerblue4", "white", "indianred3"),limits=c(-0.52,0.52),guide=guide_colorbar(frame.colour="black",ticks.colour="black",direction="horizontal"))+
    labs(x = "", y = "", fill = "Pearson\'s r") +
    theme_bw()+
    theme(legend.position = "top")
dev.off()
finalall=final[final$flag=='All Peaks',]
pdf(paste0(output_path,'AllCA_bnomialp_Peaks_chromatin_state_ASEcor_orderbydivergent.pdf'), width=3.5,height=5,useDingbats=FALSE)
ggplot(finalall,aes(factor(celltype),paper_name))+ 
    geom_tile(aes(fill=pearsonr),colour = "black")+
    geom_text(aes(label=labelest),size=2)+
    scale_fill_gradientn(colours = c("dodgerblue4", "white", "indianred3"),limits=c(-0.52,0.52),guide=guide_colorbar(frame.colour="black",ticks.colour="black",direction="horizontal"))+
    labs(x = "", y = "", fill = "Pearson\'s r") +
    theme_bw()+   
    theme(legend.position = "top")
dev.off()


