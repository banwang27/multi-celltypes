library("DESeq2")
library("qvalue")
library("stringr")
#Pairwise detection for tissue specific expressed genes
align2 <-'human'
input_path=paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',align2,'/')
# Read in raw counts for all samples
tpm<-read.csv(paste0(input_path,'rawcounts.txt'),sep='\t')
print(head(tpm))
nsamples=(ncol(tpm)-2)/3
subtpm<-tpm[,grep('total',colnames(tpm))]
output_path <-paste0('/scratch/users/banwang3/celllines/analysis/DESeq2_norm/',align2,'/')

idx <-which(rowSums(subtpm)>0)  
subtpm<-subtpm[idx,]
genenames_all<-tpm$gene[idx]
cond_line_all<-factor(str_split_fixed(colnames(subtpm), "_", 4)[,3])
sampleid<-str_split_fixed(colnames(subtpm), "_", 4)[,1]
id<-as.integer(gsub('RNA','',sampleid))
cond_batch_all<-c()
for (i in id){
    if (i<=18){
        cond_batch_all=c(cond_batch_all,1)
    }else{cond_batch_all=c(cond_batch_all,2)}
}
cond_batch_all<-factor(cond_batch_all)
cond_sample_all<- factor(id)
celltypes<-str_split_fixed(colnames(subtpm), "_", 4)[,2]
cts<-unique(celltypes)
subtpm$gene<-genenames_all
for (ct in cts){
    #keep genes that have non zero allelic counts in at least 1 sample
    ctcountsfilter<-read.csv(paste0(input_path,ct,'_totalcounts_gene_bins.txt'),sep='\t') 
    subtpm1<-subtpm[subtpm$gene%in%ctcountsfilter$gene,]
    genenames<-subtpm1$gene  
    for (ct1 in cts[cts!=ct]){
        idx=which(celltypes %in% c(ct,ct1)) 
        cond_line<-cond_line_all[idx]
        cond_batch<-cond_batch_all[idx]
        cond_sample<-cond_sample_all[idx]
        cond_ChiHum<-as.factor(celltypes[idx])
        valMat<-subtpm1[,idx]
        if (length(unique(cond_batch))==1){
            cond <- ~ cond_line + cond_ChiHum
            cond_reduced <- ~ cond_line
        }else{
            cond <- ~ cond_batch+cond_line + cond_ChiHum
            cond_reduced <- ~ cond_batch+cond_line
        }
        coldata <- data.frame(cond_batch,cond_line, cond_sample, cond_ChiHum)
        dds <- DESeqDataSetFromMatrix(valMat, coldata, cond)
        #dds <- DESeq(dds,test="LRT",betaPrior=FALSE,full=cond,reduced=~1)
        dds <- DESeq(dds,test="LRT",betaPrior=FALSE,full=cond,reduced=cond_reduced)
        res <- results(dds)
        resultsNames(dds)
        coef <- grep("cond_ChiHum_",resultsNames(dds))
        resLFC <- lfcShrink(dds, coef=coef, res=res) 
        res$padj_Mine <- p.adjust(res$pvalue, method="BH")
        resLFC$padj_Mine <- p.adjust(resLFC$pvalue, method="BH")
        #compute Storey's qvalues
        idx <- which(!is.na(resLFC$pvalue))
        resLFC$qvalue_Mine <- rep(NA,times=length(resLFC$pvalue))
        pvals_temp <- resLFC$pvalue[idx]
        qvals_temp <- qvalue(p = pvals_temp)
        resLFC$qvalue_Mine[idx] <- qvals_temp$qvalues
        resLFC$gene <- genenames
        #print outputs
        output_filename <- paste0(output_path,ct,'_vs_',ct1,'_aligned2',align2,'_qvalue_','DESeq2.txt')
        output_file <- cbind(resLFC[,"gene"],resLFC[,"log2FoldChange"],resLFC[,"pvalue"],resLFC[,"padj"],resLFC[,"padj_Mine"],resLFC[,"qvalue_Mine"])
        colnames(output_file) <- c("gene","log2FoldChange","pvalue","padj","padj_mine","qval_mine")
        write.table(output_file,file=output_filename,sep="\t", quote = FALSE,na = "", row.names = FALSE, col.names = TRUE)
    }
}
