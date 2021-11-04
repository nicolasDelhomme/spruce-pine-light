#' Working directory
setwd("/mnt/picea/projects/spruce/rgarcia/light")

#' Data
spruce.mat <- read.csv("analysis/HTSeq/spruce_raw-unnormalised_technical-replicates-merged_data.csv",
                       row.names = 1)

#' Read the annotation
annot <- read.delim("/mnt/picea/storage/reference/Picea-abies/v1.1/annotation/Pabies-Pfam-description-annotation.tsv")
colnames(annot) <- c("GeneID","Confidence","PfamDomains")
stopifnot(all(annot$GeneID == rownames(spruce.mat)))

# libraries
library(DESeq2)
library(Glimma)

#' DESeq
dds <- DESeqDataSetFromMatrix(
  countData = spruce.mat,
  colData = data.frame(
    light=factor(substr(colnames(spruce.mat),1,1)),
    ecotype=factor(substr(colnames(spruce.mat),3,5))
  ),
  design = ~ ecotype * light)

dds <- DESeq(dds)
res <- results(dds)

#' Glimma
#' ### Interactive plot
res.df <- as.data.frame(res)
res.df$log10MeanNormCount <- log10(res.df$baseMean + 1)
idx <- rowSums(counts(dds)) > 0
res.df <- res.df[idx,]
res.df$padj[is.na(res.df$padj)] <- 1
sel <- !is.na(res$padj) & res$padj <= 0.01 & abs(res$log2FoldChange) >= 0.5

#' When exploring this report, it is essential to remember that the 
#' dot plot on the left contain values that have simply been corrected
#' for their relative sequencing depth
glMDPlot(res.df,
         xval = "log10MeanNormCount",
         yval = "log2FoldChange",
         counts = counts(dds)[idx,],
         anno = annot[idx,],
         groups = dds$light,
         samples = paste(dds$ecotype,dds$light,sep="-"),
         status = sel[idx],
         display.columns = c("GeneID","Confidence","PfamDomains"),
         id.column = "GeneID",
         path = "analysis/HTSeq",
         folder = paste0("report"),
         launch=FALSE)

