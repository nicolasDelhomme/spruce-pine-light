#' ---
#' title: "Sari's Scots Pine light experiment plots"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/rgarcia/light")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/rgarcia/light")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(RColorBrewer))

#' Palettes
pal <- brewer.pal(3,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Read the data
spruce.mat <- read.csv("analysis/HTSeq/spruce_raw-unnormalised_technical-replicates-merged_data.csv",
                       row.names = 1)

pine.mat <- read.csv("analysis/HTSeq/pine_raw-unnormalised_technical-replicates-merged_data.csv",
                       row.names = 1)

#' Define an outlier
outlier <- "A.P.56N.2"

#' Rename the samples
spruceIDs <- read.delim("~/Git/UPSCb/projects/spruce-pine-light/doc/spruce-ID-mapping.txt",
                        header = FALSE, as.is = TRUE)
pineIDs <- read.delim("~/Git/UPSCb/projects/spruce-pine-light/doc/pine-ID-mapping.txt",
                      header = FALSE, as.is = TRUE)

colnames(spruce.mat) <- spruceIDs[match(colnames(spruce.mat),spruceIDs[,1]),2]
colnames(pine.mat) <- pineIDs[match(colnames(pine.mat),pineIDs[,1]),2]

#' Reorder the samples
spruce.mat <- spruce.mat[,order(colnames(spruce.mat))]
pine.mat <- pine.mat[,order(colnames(pine.mat))]

#' # Process
#' ## Variance Stabilisation
#' ### Spruce
spruce.dds <- DESeqDataSetFromMatrix(
  countData = spruce.mat,
  colData = data.frame(condition=colnames(spruce.mat)),
  design = ~ condition)

vsd <- varianceStabilizingTransformation(spruce.dds, blind=TRUE)
spruce.vst <- assay(vsd)
spruce.vst <- spruce.vst - min(spruce.vst)

write.csv(spruce.vst,"analysis/HTSeq/spruce_library-size-normalized_variance-stabilized_technical-replicates-merged_data.csv")

#' ### Pine
pine.dds <- DESeqDataSetFromMatrix(
  countData = pine.mat,
  colData = data.frame(condition=colnames(pine.mat)),
  design = ~ condition)

vsd <- varianceStabilizingTransformation(pine.dds, blind=TRUE)
pine.vst <- assay(vsd)
pine.vst <- pine.vst - min(pine.vst)

write.csv(pine.vst,"analysis/HTSeq/pine_library-size-normalized_variance-stabilized_technical-replicates-merged_data.csv")

sel <- colnames(pine.mat) != outlier
pine.dds.n.o <- DESeqDataSetFromMatrix(
  countData = pine.mat[,sel],
  colData = data.frame(condition=colnames(pine.mat[,sel])),
  design = ~ condition)

vsd <- varianceStabilizingTransformation(pine.dds.n.o, blind=TRUE)
pine.vst.n.o <- assay(vsd)
pine.vst.n.o <- pine.vst.n.o - min(pine.vst.n.o)

#' # Figures
#' ## PCA
#' ### Spruce
#' #### Including Dark
dat <- spruce.vst
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=c(17,19)[as.integer(factor(substr(colnames(dat),4,5)))],
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B","Dark","Lat. 56", "Lat. 67"),
       pch=c(15,15,15,2,1),
       col=c(pal[1:3],1,1),
       bty="n")

#' #### Excluding Dark
#' Remove the D samples (Dark)
dat <- spruce.vst[,substr(colnames(spruce.vst),1,1) != "D"]
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=c(17,19)[as.integer(factor(substr(colnames(dat),4,5)))],
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B","Lat. 56", "Lat. 67"),
       pch=c(15,15,2,1),
       col=c(pal[1:2],1,1),
       bty="n")

#' ### Pine
#' #### Including Dark
dat <- pine.vst
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=c(17,19)[as.integer(factor(substr(colnames(dat),4,5)))],
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B","Dark","Lat. 56", "Lat. 67"),
       pch=c(15,15,15,2,1),
       col=c(pal[1:3],1,1),
       bty="n")

#' #### Including Dark - no outlier
dat <- dat[,colnames(dat) != outlier]
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=c(17,19)[as.integer(factor(substr(colnames(dat),4,5)))],
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B","Dark","Lat. 56", "Lat. 67"),
       pch=c(15,15,15,2,1),
       col=c(pal[1:3],1,1),
       bty="n")

#' Not affected by the normalisation
dat <- pine.vst.n.o
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=c(17,19)[as.integer(factor(substr(colnames(dat),4,5)))],
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B","Dark","Lat. 56", "Lat. 67"),
       pch=c(15,15,15,2,1),
       col=c(pal[1:3],1,1),
       bty="n")

#' #### Excluding Dark
#' Remove the D samples (Dark)
dat <- pine.vst[,substr(colnames(pine.vst),1,1) != "D"]
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=c(17,19)[as.integer(factor(substr(colnames(dat),4,5)))],
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topleft",legend = c("A","B","Lat. 56", "Lat. 67"),
       pch=c(15,15,2,1),
       col=c(pal[1:2],1,1),
       bty="n")

#' #### Excluding Dark - no outlier
dat <- dat[,colnames(dat) != outlier]
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=c(17,19)[as.integer(factor(substr(colnames(dat),4,5)))],
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B","Lat. 56", "Lat. 67"),
       pch=c(15,15,2,1),
       col=c(pal[1:2],1,1),
       bty="n")


#' ## Heatmap
"geneSelect" <- function(cnt,splt,exp=1,nrep=2){
  rowSums(sapply(lapply(
    split.data.frame(t(cnt >= exp),splt)
    ,colSums), ">=", nrep)) >= 1
}

"repSds" <- function(vst,ord,cpu=4){
  sds <- do.call("cbind",
                 mclapply(
                   split.data.frame(t(vst),ord),colSds,
                   mc.cores=cpu))
  sds[is.na(sds)] <- 0
  return(sds)
}

#' ### Spruce
#' Select the genes that have vst >= 4 in at least 2 rep out of 3 where the sd within rep is <= 1
sel <- geneSelect(spruce.vst,
                  sub("\\.\\d+$","",colnames(spruce.vst)),
                  exp=4)

sds <- repSds(spruce.vst,sub("\\.\\d+$","",colnames(spruce.vst)))
sel <- sel & rowSums(sds <= 1) == 6
message(sprintf("There are %s surviving the cutoffs",sum(sel)))

#' #### Including Dark
#' expression values
dat <- spruce.vst[sel,]
heatmap.2(as.matrix(dat),
          labRow=NA,trace="none",
          las=2,col=hpal)

#' saturated expression
dat.sat <- dat
dat.sat[dat < 3 ] <- 3
dat.sat[dat > 8 ] <- 8
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal)

#' standard scores
heatmap.2(as.matrix(dat),
          scale="row",
          labRow=NA,trace="none",
          las=2,col=hpal)

#' #### Excluding Dark
dat <- dat[,substr(colnames(dat),1,1) != "D"]
#' expression values
heatmap.2(as.matrix(dat),
          labRow=NA,trace="none",
          las=2,col=hpal,cexCol = 0.8)

#' saturated expression
dat.sat <- dat
dat.sat[dat < 3 ] <- 3
dat.sat[dat > 8 ] <- 8
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,cexCol = 0.8)

#' saturated expression, no sample clustering
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,Colv = NULL,
          dendrogram="row",cexCol = 0.8)

#' standard scores
heatmap.2(as.matrix(dat),
          scale="row",
          labRow=NA,trace="none",
          las=2,col=hpal,cexCol = 0.8)

#' ### Pine
#' Select the genes that have vst >= 2 in at least 2 rep out of 3 where the sd within rep is <= 1 in at least
#' 5 out of the 6 set of replicates

#' Using the outlier or not has little effect
sel.n.o <- geneSelect(pine.vst.n.o,
                      sub("\\.\\d+$","",colnames(pine.vst.n.o)),
                      exp=2)
sds <- repSds(pine.vst.n.o,sub("\\.\\d+$","",colnames(pine.vst.n.o)))
plot(density(sds))
sel <- sel & rowSums(sds <= 1) >= 5
message(sprintf("There are %s surviving the cutoffs",sum(sel)))

#' Including the outlier
sel <- geneSelect(pine.vst,
                  sub("\\.\\d+$","",colnames(pine.vst)),
                  exp=2)

sds <- repSds(pine.vst,sub("\\.\\d+$","",colnames(pine.vst)))
plot(density(sds))
sel <- sel & rowSums(sds <= 1) >= 5
message(sprintf("There are %s surviving the cutoffs",sum(sel)))

#' #### Including Dark
#' expression values
dat <- pine.vst[sel,]
heatmap.2(as.matrix(dat),
          labRow=NA,trace="none",
          las=2,col=hpal)

#' saturated expression
dat.sat <- dat
dat.sat[dat < 2 ] <- 2
dat.sat[dat > 8 ] <- 8
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal)

#' standard scores
heatmap.2(as.matrix(dat),
          scale="row",
          labRow=NA,trace="none",
          las=2,col=hpal)

#' #### Excluding Dark
#' expression values
dat <- dat[,substr(colnames(dat),1,1) != "D"]
heatmap.2(as.matrix(dat),
labRow=NA,trace="none",
las=2,col=hpal,cexCol = 0.8)

#' saturated expression
dat.sat <- dat
dat.sat[dat < 2 ] <- 2
dat.sat[dat > 8 ] <- 8
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,cexCol = 0.8)

#' saturated expression, no sample clustering
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,Colv = NULL,
          dendrogram="row",cexCol = 0.8)

#' standard scores
heatmap.2(as.matrix(dat),
          scale="row",
          labRow=NA,trace="none",
          las=2,col=hpal,cexCol = 0.8)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
