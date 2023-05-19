
library(DESeq2)
library(GenomicRanges)
library(qvalue)
library(apeglm)
library(dplyr)
library(ggplot2)
library(gplots)
library(gridExtra)
#library(Rtsne)
library(reshape)
library(scales)
#library(VennDiagram)
#library(Seurat)
library(graphics)
#library(MultiPhen)
library(stringr)

#if (length(args)<6) {
#  stop("Need to supply config file, ase (boolean), indir (directory for files)", call.=FALSE)
#}
#Read in the config file and the ASE file
config_DESeq <- 0
for (x in list.files()) {
  if (grepl("Config", x, fixed = TRUE) & grepl("27", x, fixed = TRUE)) {
    config_DESeq <- x
    print(config_DESeq)
    ase <- FALSE
    indir <- "C:/Users/astar/PrimeDB/"
    outfile = str_replace(config_DESeq, "_Config", "_DESeq2")
    
    config <- read.csv(config_DESeq, sep = '\t', header = FALSE)
    files <- config$V1
    cond_of_interest <- config$V2
    l <- length(config)
    
    #Infer the two types being compared.
    #Assumes that human comes first.
    if (ase == FALSE) {
      type_1 <- unique(cond_of_interest)[1]
      type_2 <- unique(cond_of_interest)[2]
    }
    
    #Read in the genes and only keep those present in both.
    valMat = c()
    valMat_1 = c()
    valMat_2 = c()
    for (rep in 1:length(files)) {
      input_filename <- files[rep]
      repVals <- read.table(input_filename , header=F, na.strings = "", sep='\t', stringsAsFactors=F)
      if (cond_of_interest[rep] == type_1) {
        x1 <- repVals$V1
      } else {
        x2 <- repVals$V1
      }
    }
    
    keep <- intersect(x1, x2)
    keep <- keep[! keep %in% c('__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique')]

    
    #For each file, compute the TPM. 
    a <- c()
    t <- c()
    #Alter this to only compute for genes of interest (in z)
    for (file in files) {
      f <- read.table(file, header = FALSE, sep = "\t")
      rownames(f) = f$V1
      f <- subset(f, f$V1 %in% keep)
      f <- dplyr::select(f,-V1)
      if (grepl("Human", file, fixed=TRUE))  or (grepl("human", file, fixed=TRUE)){
        geneinfo = "human_gene_length.txt"
      } else if (grepl("Chimp", file, fixed=TRUE))  or (grepl("chimp", file, fixed=TRUE)){
        geneinfo = "chimp_gene_length.txt"
      } else if (grepl("Rhesus", file, fixed=TRUE)) or (grepl("rhesus", file, fixed=TRUE)) {
        geneinfo = "rhesus_gene_length.txt"
      } else if (grepl("Orangutan", file, fixed=TRUE)) or (grepl("orangutan", file, fixed=TRUE)) {
        geneinfo = "orangutan_gene_length.txt"
      } else if (grepl("Mouse", file, fixed=TRUE)) or (grepl("mouse", file, fixed=TRUE)) {
        geneinfo = "mouse_gene_length.txt"
      } else if (grepl("Gorilla", file, fixed=TRUE)) or (grepl("gorilla", file, fixed=TRUE)) {
        geneinfo = "gorilla_gene_length.txt"
      } else if (grepl("Marmoset", file, fixed=TRUE)) or (grepl("marmoset", file, fixed=TRUE)) {
        geneinfo = "marmoset_gene_length.txt"
      } else if (grepl("Bonobo", file, fixed=TRUE)) or (grepl("bonobo", file, fixed=TRUE)) {
        geneinfo = "bonobo_gene_length.txt"
      }
      geneInfo = read.table(geneinfo, header=TRUE)
      geneInfo$genesizeKb = as.vector(as.numeric(geneInfo$merged)/1000)
      f_tot_cpm = sweep(as.matrix(f), 2, as.double(colSums(f)/1000000), `/`)
      f_tot_fpkm = sweep(as.matrix(f_tot_cpm), 1, as.double(geneInfo$genesizeKb), `/`)
      f_tot_rpk = sweep(as.matrix(f), 1, as.double(geneInfo$genesizeKb), `/`)
      f_tot_tpm = sweep(as.matrix(f_tot_rpk), 2, as.double(colSums(f_tot_rpk)/1000000), `/`)
      a <- cbind(a, as.data.frame(f_tot_tpm)$V2)
      t <- c(t, file)
    }
    rownames(a) <- rownames(f_tot_cpm)
    colnames(a) <- t
    write.table(file = paste(substr(outfile,1,nchar(file)-4), "_tpm.txt", sep = ""),
                a, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    #From this github page: https://github.com/konopkalab/Primates_CellType/blob/master/DGE/Dge_NeuN_LM.R
    ## Variance Explained
    # counts = gene x sample matrix
    # meta = sample x covariate matrix
    # threshold = number of PCA to consider (e.g. 5)
    # inter = interaction between covariates (e.g. Age:Sex)
    VarExp <- function(counts, meta, threshold, inter){
      suppressPackageStartupMessages(library(lme4))
      suppressPackageStartupMessages(library(optimx))
      counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
      cor.counts <- cor(counts.center)
      dim(cor.counts)
      eigen.counts <- eigen(cor.counts)
      eigen.mat <- eigen.counts$vectors
      eigen.val <- eigen.counts$values
      n.eigen <- length(eigen.val)
      eigen.val.sum <- sum(eigen.val)
      percents.pcs <- eigen.val/eigen.val.sum
      meta <- as.data.frame(meta)
      
      all <- 0
      npc.in <- 0
      for(i in 1:n.eigen){
        all <- all + percents.pcs[i]
        npc.in <- npc.in + 1
        if(all > threshold){break}
      }
      if (npc.in < 3) {npc <- 3}
      
      pred.list <- colnames(meta)
      meta <- droplevels(meta)
      
      n.preds <- ncol(meta) + 1
      if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}
      
      ran.pred.list <- c()
      for(i in 1:ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
      }
      ##interactions
      if(inter){
        for(i in 1:(ncol(meta)-1)){
          for(j in (i+1):ncol(meta)){
            ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
            pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
          }
        }
      }
      formula <- paste(ran.pred.list, collapse = " + ")
      formula <- paste("pc", formula, sep=" ~ ")
      ran.var.mat <- NULL
      for(i in 1:npc.in){
        dat <- cbind(eigen.mat[,i],meta)
        colnames(dat) <- c("pc",colnames(meta))
        Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit,control = lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
        var.vec <- unlist(VarCorr(Rm1ML))
        ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
      }
      ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
      wgt.vec <- eigen.val/eigen.val.sum
      prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
      std.prop.var <- prop.var/sum(prop.var)
      std.prop.var
    }
    
    # Bar plot for data visualization
    plotVarExp <- function(pvca.res, title){
      suppressPackageStartupMessages(library(ggplot2))
      plot.dat <- data.frame(eff=names(pvca.res), prop=pvca.res)
      p <- ggplot(plot.dat, aes(x=eff, y=prop))+
        ggtitle(title)+
        geom_bar(stat="identity", fill="steelblue", colour="steelblue") +
        geom_text(aes(label=round(prop,3), y=prop+0.04), size=4) +
        scale_x_discrete(limits=names(pvca.res)) +
        scale_y_continuous(limits = c(0,1)) +
        labs(x= "Effects", y= "WAPV") +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      p
    }
    
    #Will use uncorrected TPM a la GTEx, but will compute the variation explained by each covariate yay?
    if (l == 3) {
      cov1 = config$V3
      pdf(paste("Variance_Explained_", outfile, ".pdf", sep = ""),width=8,height=6)
      bleh <- cbind(cond_of_interest, cov1)
      print(bleh)
      var <- VarExp(a, bleh, l-1, FALSE)
      print(plotVarExp(var,"Variance Explained"))
      dev.off()
    } else if (l == 4) {
      cov1 = config$V3
      cov2 = config$V4
      pdf(paste("Variance_Explained_", outfile, ".pdf", sep = ""),width=8,height=6)
      bleh <- cbind(cond_of_interest, cov1, cov2)
      print(bleh)
      var <- VarExp(a, bleh, l-1, FALSE)
      print(plotVarExp(var,"Variance Explained"))
      dev.off()
    } else if (l == 5) {
      cov1 = config$V3
      cov2 = config$V4
      cov3 = config$V5
      pdf(paste("Variance_Explained_", outfile, ".pdf", sep = ""),width=8,height=6)
      bleh <- cbind(cond_of_interest, cov1, cov2, cov3)
      print(bleh)
      var <- VarExp(a, bleh, l-1, FALSE)
      print(plotVarExp(var,"Variance Explained"))
      dev.off()
    } else if (l == 6) {
      cov1 = config$V3
      cov2 = config$V4
      cov3 = config$V5
      cov4 = config$V6
      pdf(paste("Variance_Explained_", outfile, ".pdf", sep = ""),width=8,height=6)
      bleh <- cbind(cond_of_interest, cov1, cov2, cov3, cov4)
      print(bleh)
      var <- VarExp(a, bleh, l-1, FALSE)
      print(plotVarExp(var,"Variance Explained"))
      dev.off()
    } else if (l == 7) {
      cov1 = config$V3
      cov2 = config$V4
      cov3 = config$V5
      cov4 = config$V6
      cov5 = config$V7
      pdf(paste("Variance_Explained_", outfile, ".pdf", sep = ""),width=8,height=6)
      bleh <- cbind(cond_of_interest, cov1, cov2, cov3, cov4, cov5)
      print(bleh)
      var <- VarExp(a, bleh, l-1, FALSE)
      print(plotVarExp(var,"Variance Explained"))
      dev.off()
    } else if (l > 7) {
      print("Too Many Covariates, I was lazy and hardcoded up to 5")
    } else {
      print("No Covariates")
      pdf(paste("Variance_Explained_", outfile, ".pdf", sep = ""),width=8,height=6)
      bleh <- cbind(cond_of_interest)
      print(bleh)
      var <- VarExp(a, bleh,l-1, FALSE)
      print(plotVarExp(var,"Variance Explained"))
      dev.off()
    }
    
    
  }
}
