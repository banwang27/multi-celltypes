library(DESeq2)
library(stringr)
library(qvalue)
Aligned2 <-'human'
input_path <- paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',Aligned2,'/')
endFile <- '_ase_normalized.txt'
output_path_normalized <-paste0('/scratch/users/banwang3/celllines/analysis/ASE_norm/',Aligned2,'/deseq2_normalized_counts/')
output_path <-paste0('/scratch/users/banwang3/celllines/analysis/DESeq2_norm/',Aligned2,'/')
all_files = list.files(input_path)
sub_all_files = all_files[grep(endFile, all_files)]
celltypes=unique(str_split_fixed(sub_all_files, "_", 4)[,2])
for (cell_type in celltypes){
    sub_files=sub_all_files[grep(cell_type,sub_all_files)]
    valMat = c()
    #allhumangenes<-read.table('/scratch/users/banwang3/celllines/analysis/ASE/human/withPseudogenes/All_Human_Genes.csv',sep=',',header=TRUE)
    for (a in sub_files){  
        input_filename<-paste0(input_path,a)
        repVals <- read.table(input_filename , header=TRUE, na.strings = "", sep='\t', stringsAsFactors=F)
        valMat <- cbind(valMat,repVals[,3])
        valMat <- cbind(valMat,repVals[,4])
    }
    genenames <-repVals$gene
    #Filter out genes have no allele expression
    idx <-which(rowSums(valMat)>0)  
    valMat<-valMat[idx,]
    genenames<-genenames[idx]
    num_samples = ncol(valMat) / 2
    cond_sample <- factor(rep(1:num_samples,each=2)) 
    cond_ChiHum <- factor(rep(c("Human","Chimp"),times=num_samples)) #first level
    cond_line<-factor(rep(str_split_fixed(sub_files, "_", 4)[,3],each=2))
    sampleid<-str_split_fixed(sub_files, "_", 4)[,1]
    #Find out batch effect, for sample have ID greater than 18 is from batch 2
    id<-as.integer(gsub('RNA','',sampleid))
    cond_batch<-c()
    for (i in id){
        if (i<=18){
            cond_batch=c(cond_batch,1,1)
        }else{cond_batch=c(cond_batch,2,2)}
    }
    cond_batch<-factor(cond_batch)
    #Find out if there is batch effect in this cell type data
    if (length(unique(cond_batch))==1){
        cond <- ~ cond_line + cond_ChiHum
        cond_reduced <- ~ cond_line
    }else{
        cond <- ~ cond_batch+cond_line + cond_ChiHum
        cond_reduced <- ~ cond_batch+cond_line
    }
    coldata <- data.frame(cond_batch,cond_line, cond_sample, cond_ChiHum)
    #Construct DESeq object and normalize read counts (median of ratios)
    dds <- DESeqDataSetFromMatrix(valMat, coldata, cond)
    dds <- estimateSizeFactors(dds)
    normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
    colnames(normalized_counts)=paste0(rep(gsub('_ase_normalized.txt','',sub_files),each=2),'_',cond_ChiHum)
    normalized_counts$gene=genenames
    # Save the normalized read counts
    write.table(normalized_counts, file=paste0(output_path_normalized,cell_type,'_','deseq2_normalized_counts.txt'), sep="\t", quote=F, col.names=TRUE,row.names = FALSE)
    # Find DE genes using DESeq2
    dds <- DESeq(dds,test="LRT",betaPrior=FALSE,full=cond,reduced=cond_reduced)
    res <- results(dds)
    print(resultsNames(dds))
    coef <- grep("Human_vs_",resultsNames(dds))
    resLFC <- lfcShrink(dds, coef=coef, res=res) 
    res$padj_Mine <- p.adjust(res$pvalue, method="BH")
    resLFC$padj_Mine <- p.adjust(resLFC$pvalue, method="BH")
    #compute Storey's qvalues
    idx <- which(!is.na(resLFC$pvalue))
    resLFC$qvalue_Mine <- rep(NA,times=length(resLFC$pvalue))
    pvals_temp <- resLFC$pvalue[idx]
    qvals_temp <- qvalue(p = pvals_temp)
    resLFC$qvalue_Mine[idx] <- qvals_temp$qvalues
    filename <- paste0(output_path,cell_type,'_aligned2',Aligned2,'_MA.pdf')
    pdf(filename, useDingbats=FALSE)
    print(plotMA(res))
    dev.off()
    filename <- paste0(output_path,cell_type,'_aligned2',Aligned2,'_MA_shrunken.pdf')
    pdf(filename, useDingbats=FALSE)
    print(plotMA(resLFC))
    dev.off()
    
    resLFC$gene <- genenames
    #print outputs
    output_filename <- paste0(output_path,cell_type,'_aligned2',Aligned2,'_qvalue_','DESeq2.txt')
    output_file <- cbind(resLFC[,"gene"],resLFC[,"log2FoldChange"],resLFC[,"pvalue"],resLFC[,"padj"],resLFC[,"padj_Mine"],resLFC[,"qvalue_Mine"])
    colnames(output_file) <- c("gene","log2FoldChange","pvalue","padj","padj_mine","qval_mine")
    write.table(output_file,file=output_filename,sep="\t", quote = FALSE,na = "", row.names = FALSE, col.names = TRUE)
}
