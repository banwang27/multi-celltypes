library(stringr)
library(ggplot2)

getCorString <- function(dataframe, x, y, category_col = NULL, labelRound = 2, method = "pearson") {
    # Calculates pearson correlation coefficient and pvalue between two columns in a dataframe
    # Construct category list of factors
    if(!is.null(category_col)) {
        catList <- list()
        for(i in category_col){
            catList[[i]] <- dataframe[,i]
           }
    } else {
         catList <- rep(T, nrow(dataframe))
    }
    # Loop by categorie
    results <- by(dataframe, catList, function(x, xVal, yVal, category_col) { 
                  # Get correlations
                  corR <- cor.test(x[,xVal], x[,yVal], method = method)
                  corText <- paste0("r = ", signif(corR$estimate,labelRound) , "\np = ", signif(corR$p.value, labelRound))
                  result <- data.frame(r = corR$estimate, p =  corR$p.value, text = corText, stringsAsFactors = F)
                  # Append category names
                  if(!is.null(category_col)) {
                      for(i in category_col){
                       result <- cbind (result, x[1,i], stringsAsFactors = F)
                      }
                      colnames(result)[ (ncol(result) - length(category_col) + 1) : ncol(result) ] <- category_col
                   }
                   return(result)
                }, xVal = x, yVal = y, category_col = category_col)
     results <- do.call(rbind, results)
     return(results)
}
input_path='/scratch/users/banwang3/celllines/ATAC/findSNPs/summarytable/'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/suppplots/'

df = read.csv(paste0(input_path,'fasta_sei_pred_all34SNPsin4096bp_lfc.txt'),sep='\t')
idx = which(df[,paste0('humanfasta')]>=0.5| df[,paste0('chimpfasta')]>=0.5)
df1 = df[idx,]

human=c()
for (i in seq(1,34)){ 
    sei = data.frame(lfc = df1[,paste0('humanbased_humchplfc_SNP',i)],SNP = paste0('SNP',i),type='humanbasedSNP')
    human=rbind(human,sei)
}
sei = data.frame(lfc = df1[,'fastalfc'],SNP='Fasta',type='fasta')
human = rbind(human,sei)
human$SNP = factor(human$SNP,levels=c('Fasta',paste0('SNP',seq(1,34))))

pdf(paste0(output_path,'boxplot_fasta_humanbased_lfc_all34SNPs_centeratHAR.pdf'),useDingbats=FALSE,width=15,height=4)
ggplot(human,aes(x = SNP,y=lfc,fill=type))+
     geom_hline(yintercept=0, color = "grey")+
    geom_boxplot(outlier.size = 0.15,show.legend=FALSE)+
     ylab('Log2(Human/Chimp)')+
     xlab('')+
     theme_bw()+
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()

data = data.frame(x=human[human$SNP=='Fasta','lfc'],y=human[human$SNP=='SNP21','lfc'])
texts=getCorString(data,'x','y')
texts$x=-Inf
texts$y=Inf
pdf(paste0(output_path,'scatterplot_fasta_humanbased_lfc_centeratHAR_SNP21.pdf'),useDingbats=FALSE,width=3,height=3)
ggplot(data,aes(x=x,y=y))+
    geom_smooth(method=lm,se=FALSE, fullrange=TRUE,linewidth=0.5,alpha=0.2,linetype = "longdash")+
    geom_point(size=0.2,alpha=0.3,show.legend = FALSE)+
    geom_text(data=texts,aes(label=text),hjust=0,vjust=1,size=3)+
    xlab('Log2(Human/Chimp) - Fasta')+
    ylab('Log2(Human/Chimp) - SNP21')+
    theme_bw()+
    theme(axis.title = element_text(size = 8),legend.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
 
