#' ---
#' title: "Sari's Scots Pine light experiment PCA and Heatmap"
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
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(RColorBrewer))

#' Helper
source("~/Git/UPSCb/src/R/featureSelection.R")
source("~/Git/UPSCb/src/R/expressionSpecificityUtility.R")

#' Palettes
pal <- brewer.pal(3,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Create the list of genes of interest
gene.list <- mclapply(dir("~/Git/UPSCb/projects/spruce-pine-light/doc",
             pattern="-1percent.csv$",
             full.names = TRUE),
         function(f){
           co <- read.csv(f,as.is = TRUE)[,c("ID_pine","ID_spruce")]
           list(pine=unique(co$ID_pine[co$ID_pine != ""]),
                spruce=unique(co$ID_spruce[co$ID_spruce != ""]))
         },mc.cores=4L)

p.genes <- unique(unlist(sapply(gene.list,"[","pine"),use.names = FALSE))
p.genes <- sub("\\.0$","",p.genes)
s.genes <- unique(unlist(sapply(gene.list,"[","spruce"),use.names = FALSE))
  
#' Read the data
spruce.mat <- read.csv("analysis/HTSeq/spruce_raw-unnormalised_technical-replicates-merged_data.csv",
                       row.names = 1)

pine.mat <- read.csv("analysis/HTSeq/pine_raw-unnormalised_technical-replicates-merged_data.csv",
                       row.names = 1)
rownames(pine.mat) <- sub("\\.0$","",rownames(pine.mat))

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

#' Select only the samples from 67N
spruce.mat <- spruce.mat[,grep("67N",colnames(spruce.mat))]
pine.mat <- pine.mat[,grep("67N",colnames(pine.mat))]

#' # Process
#' ## Variance Stabilisation - blind
#' ### Spruce
spruce.dds <- DESeqDataSetFromMatrix(
  countData = spruce.mat,
  colData = data.frame(condition=colnames(spruce.mat)),
  design = ~ condition)

vsd <- varianceStabilizingTransformation(spruce.dds, blind=TRUE)
spruce.vst <- assay(vsd)
spruce.vst <- spruce.vst - min(spruce.vst)

write.csv(spruce.vst,"analysis/HTSeq/spruce_library-size-normalized_variance-stabilized_technical-replicates-merged_67N_data.csv")

#' ### Pine
pine.dds <- DESeqDataSetFromMatrix(
  countData = pine.mat,
  colData = data.frame(condition=colnames(pine.mat)),
  design = ~ condition)

vsd <- varianceStabilizingTransformation(pine.dds, blind=TRUE)
pine.vst <- assay(vsd)
pine.vst <- pine.vst - min(pine.vst)

write.csv(pine.vst,"analysis/HTSeq/pine_library-size-normalized_variance-stabilized_technical-replicates-merged_67N_data.csv")

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
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B","Dark"),
       pch=15,
       col=pal[1:3],
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
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B"),
       pch=15,
       col=pal[1:2],
       bty="n")

#' #### Significant genes only
dat <- spruce.vst[s.genes,
                  substr(colnames(spruce.vst),1,1) != "D"]
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B"),
       pch=15,
       col=pal[1:2],
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
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B","Dark"),
       pch=15,
       col=pal[1:3],
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
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B"),
       pch=15,
       col=pal[1:2],
       bty="n")

#' #### Significant genes only
dat <- pine.vst[p.genes,
                  substr(colnames(pine.vst),1,1) != "D"]
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B"),
       pch=15,
       col=pal[1:2],
       bty="n")

#' ## Heatmap
#' ### Feature selection - Spruce
#' We select features to reduce the biological noise for visualisation
s.sels <- sapply(1:10, function(i){
  featureSelect(spruce.vst,
                factor(sub("\\.\\d+$","",colnames(spruce.vst))),
                i)
})

plot(colSums(s.sels),main="gene selected at a threshold",xlab="threshold",
     ylab="number of genes")
abline(v=4,lty=2,col="grey")

#' Select the genes that have vst >= 4 in at least 2 rep out of 3
s.sel <- s.sels[,4]

#' ### Feature selection - Pine
#' We select features to reduce the biological noise for visualisation
p.sels <- sapply(1:10, function(i){
  featureSelect(pine.vst,
                factor(sub("\\.\\d+$","",colnames(pine.vst))),
                i)
})

plot(colSums(p.sels),main="gene selected at a threshold",xlab="threshold",
     ylab="number of genes")
abline(v=2,lty=2,col="grey")


#' Select the genes that have vst >= 3 in at least 2 rep out of 3
p.sel <- p.sels[,2]

#' ### Spruce
#' #### Including Dark
#' Standard score, pearson distance, Ward clustering
s.dat <- t(scale(t(spruce.vst[s.sel,])))
heatmap.2(s.dat,labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(s.dat)))

#' expression values
dat <- spruce.vst[s.sel,]
heatmap.2(as.matrix(dat),
          labRow=NA,trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(dat)))

#' saturated expression
dat.sat <- dat
dat.sat[dat < 3 ] <- 3
dat.sat[dat > 8 ] <- 8
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(dat)))

#' #### Excluding Dark
#' Standard score, pearson distance, Ward clustering
s.dat <- s.dat[,substr(colnames(s.dat),1,1) != "D"]
heatmap.2(s.dat,labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(s.dat)))

#' expression values
dat <- dat[,substr(colnames(dat),1,1) != "D"]
heatmap.2(as.matrix(dat),
          labRow=NA,trace="none",
          las=2,col=hpal,cexCol = 0.8,
          labCol = sub("\\.S\\.67N\\.","_",colnames(dat)))

#' saturated expression
dat.sat <- dat
dat.sat[dat.sat < 3 ] <- 3
dat.sat[dat.sat > 8 ] <- 8
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,cexCol = 0.8,
          labCol = sub("\\.S\\.67N\\.","_",colnames(dat)))

#' saturated expression, no sample clustering
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,Colv = NULL,
          dendrogram="row",cexCol = 0.8,
          labCol = sub("\\.S\\.67N\\.","_",colnames(dat)))

#' #### Gene of interest
s.dat <- t(scale(t(spruce.vst[s.genes,
                              substr(colnames(s.dat),1,1) != "D"])))

message(sprintf("There are %s genes of interest that are not expressed",
                sum(is.na(rowSums(s.dat)))))

heatmap.2(s.dat[!is.na(rowSums(s.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(s.dat)))

#' saturated expression
dat.sat <- spruce.vst[s.genes,
                      substr(colnames(spruce.vst),1,1) != "D"]
dat.sat[dat.sat < 3 ] <- 3
dat.sat[dat.sat > 8 ] <- 8
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(s.dat)))

#' ### Pine
#' #### Including Dark
#' Standard score, pearson distance, Ward clustering
s.dat <- t(scale(t(pine.vst[p.sel,])))
heatmap.2(s.dat,labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(s.dat)))

#' expression values
dat <- pine.vst[p.sel,]
heatmap.2(as.matrix(dat),
          labRow=NA,trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(dat)))

#' saturated expression
dat.sat <- dat
dat.sat[dat < 2 ] <- 2
dat.sat[dat > 8 ] <- 8
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(dat)))

#' #### Excluding Dark
#' remove non expressed genes
dat <- dat[,substr(colnames(dat),1,1) != "D"]
dat <-dat[rowSums(dat)!=0,]

#' Standard score, pearson distance, Ward clustering
s.dat <- t(scale(t(dat)))
heatmap.2(s.dat,labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D2")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(s.dat)))

#' expression values
heatmap.2(as.matrix(dat),
labRow=NA,trace="none",
las=2,col=hpal,cexCol = 0.8,
labCol = sub("\\.P\\.67N\\.","_",colnames(dat)))

#' saturated expression
dat.sat <- dat
dat.sat[dat < 2 ] <- 2
dat.sat[dat > 8 ] <- 8
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,cexCol = 0.8,
          labCol = sub("\\.P\\.67N\\.","_",colnames(dat)))

#' saturated expression, no sample clustering
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,Colv = NULL,
          dendrogram="row",cexCol = 0.8,
          labCol = sub("\\.P\\.67N\\.","_",colnames(dat)))

#' #### Gene of interest
s.dat <- t(scale(t(pine.vst[p.genes,
                              substr(colnames(s.dat),1,1) != "D"])))

message(sprintf("There are %s genes of interest that are not expressed",
                sum(is.na(rowSums(s.dat)))))

heatmap.2(s.dat[!is.na(rowSums(s.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(s.dat)))

#' ## Variance Stabilisation - model aware
#' ### Spruce
spruce.dds <- DESeqDataSetFromMatrix(
  countData = spruce.mat,
  colData = data.frame(light=substr(colnames(spruce.mat),1,1)),
  design = ~ light)

vsd <- varianceStabilizingTransformation(spruce.dds, blind=FALSE)
spruce.vst <- assay(vsd)
spruce.vst <- spruce.vst - min(spruce.vst)

write.csv(spruce.vst,"analysis/HTSeq/spruce_library-size-normalized_variance-stabilized_technical-replicates-merged_model-aware_67N_data.csv")

#' ### Pine
pine.dds <- DESeqDataSetFromMatrix(
  countData = pine.mat,
  colData = data.frame(light=substr(colnames(pine.mat),1,1)),
  design = ~ light)

vsd <- varianceStabilizingTransformation(pine.dds, blind=FALSE)
pine.vst <- assay(vsd)
pine.vst <- pine.vst - min(pine.vst)

write.csv(pine.vst,"analysis/HTSeq/pine_library-size-normalized_variance-stabilized_technical-replicates-merged_model-aware_67N_data.csv")

#' ## Expression Specificity
#' ### Spruce
#' #### Including Dark
sES <- expressionSpecificity(spruce.vst,
                      tissues = substr(colnames(spruce.vst),1,1),
                      output = "complete")

#' #### Excluding Dark
sel <- substr(colnames(spruce.vst),1,1) != "D"
sESnD <- expressionSpecificity(spruce.vst[,sel],
                             tissues = substr(colnames(spruce.vst),1,1)[sel],
                             output = "complete")

#' ### Pine
#' #### Including Dark
pES <- expressionSpecificity(pine.vst,
                             tissues = substr(colnames(pine.vst),1,1),
                             output = "complete")

#' #### Excluding Dark
sel <- substr(colnames(pine.vst),1,1) != "D"
pESnD <- expressionSpecificity(pine.vst[,sel],
                               tissues = substr(colnames(pine.vst),1,1)[sel],
                               output = "complete")

#' ## PCA
#' ### Spruce
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
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B"),
       pch=15,
       col=pal[1:2],
       bty="n")

#' #### Significant genes only
dat <- spruce.vst[s.genes,
                  substr(colnames(spruce.vst),1,1) != "D"]
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B"),
       pch=15,
       col=pal[1:2],
       bty="n")

#' ### Pine
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
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B"),
       pch=15,
       col=pal[1:2],
       bty="n")

#' #### Significant genes only
dat <- pine.vst[p.genes,
                  substr(colnames(pine.vst),1,1) != "D"]
pc <- prcomp(t(dat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(dat),1,1)))],
     pch=17,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",legend = c("A","B"),
       pch=15,
       col=pal[1:2],
       bty="n")

#' ## Heatmap
#' ### Feature selection
s.sels <- sapply(1:10, function(i){
  featureSelect(spruce.vst,
                factor(sub("\\.\\d+$","",colnames(spruce.vst))),
                i)
})

plot(colSums(s.sels),main="gene selected at a threshold",xlab="threshold",
     ylab="number of genes")
abline(v=2,lty=2,col="grey")

#' Select the genes that have vst >= 2 in at least 2 rep out of 3
s.sel <- s.sels[,2]

#' ### Feature selection - Pine
#' We select features to reduce the biological noise for visualisation
p.sels <- sapply(1:10, function(i){
  featureSelect(pine.vst,
                factor(sub("\\.\\d+$","",colnames(pine.vst))),
                i)
})


plot(colSums(p.sels),main="gene selected at a threshold",xlab="threshold",
     ylab="number of genes")
abline(v=1,lty=2,col="grey")

#' Select the genes that have vst >= 1 in at least 2 rep out of 3
p.sel <- p.sels[,1]

#' ### Spruce
#' Standard score, pearson distance, Ward clustering
s.dat <- t(scale(t(spruce.vst[s.sel,])))
heatmap.2(s.dat,labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark)
s.dat <- t(scale(t(spruce.vst[s.genes,])))

message(sprintf("There are %s genes of interest that are not expressed",
                sum(is.na(rowSums(s.dat)))))

heatmap.2(s.dat[!is.na(rowSums(s.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark, sorted by tissue spec)
exp.peak <- c("A","B","D")[sapply(
  apply(sES[rownames(s.dat),c("A","B","D")] == sES[rownames(s.dat),"maxn"],1,which),"[",1)]
ord <- order(exp.peak,(1-sES[rownames(s.dat),"score"]))

heatmap.2(s.dat[ord,][!is.na(rowSums(s.dat[ord,])),],
          labRow=NA, dRowv = NULL, dendrogram = "column",
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark, mean of replicates)
mean.dat <- sapply(split.data.frame(t(s.dat),
                                    substr(colnames(spruce.vst),1,1)),colSums)

heatmap.2(mean.dat[!is.na(rowSums(mean.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal)

#' #### Gene of interest (including Dark, sorted by tissue spec, 
#' mean of replicates)
heatmap.2(mean.dat[ord,][!is.na(rowSums(mean.dat[ord,])),],
          labRow=NA, dRowv = NULL, dendrogram = "column",
          trace="none",
          las=2,col=hpal)

#' #### Gene of interest (excluding Dark)
nd.sel <- substr(colnames(spruce.vst),1,1)!="D"
s.dat <- s.dat[,nd.sel]

message(sprintf("There are %s genes of interest that are not expressed",
                sum(is.na(rowSums(s.dat)))))

heatmap.2(s.dat[!is.na(rowSums(s.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark, sorted by tissue spec)
exp.peak <- c("A","B")[sapply(
  apply(sESnD[rownames(s.dat),c("A","B")] == sESnD[rownames(s.dat),"maxn"],1,which),"[",1)]
ord <- order(exp.peak,(1-sESnD[rownames(s.dat),"score"]))

heatmap.2(s.dat[ord,][!is.na(rowSums(s.dat[ord,])),],
          labRow=NA, dRowv = NULL, dendrogram = "column",
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.S\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark, mean of replicates)
mean.dat <- mean.dat[,colnames(mean.dat)!="D"]

heatmap.2(mean.dat[!is.na(rowSums(mean.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal)

#' #### Gene of interest (including Dark, sorted by tissue spec, 
#' mean of replicates)
heatmap.2(mean.dat[ord,][!is.na(rowSums(mean.dat[ord,])),],
          labRow=NA, dRowv = NULL, dendrogram = "column",
          trace="none",
          las=2,col=hpal)

#' ### Pine
#' Standard score, pearson distance, Ward clustering
s.dat <- t(scale(t(pine.vst[p.sel,])))
heatmap.2(s.dat,labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark)
s.dat <- t(scale(t(pine.vst[p.genes,])))

message(sprintf("There are %s genes of interest that are not expressed",
                sum(is.na(rowSums(s.dat)))))

heatmap.2(s.dat[!is.na(rowSums(s.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark, sorted by tissue spec)
exp.peak <- c("A","B","D")[sapply(
  apply(pES[rownames(s.dat),c("A","B","D")] == pES[rownames(s.dat),"maxn"],
        1,which),"[",1)]
ord <- order(exp.peak,(1-pES[rownames(s.dat),"score"]))

heatmap.2(s.dat[ord,][!is.na(rowSums(s.dat[ord,])),],
          labRow=NA, dRowv = NULL, dendrogram = "column",
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark, mean of replicates)
mean.dat <- sapply(split.data.frame(t(s.dat),
                                    substr(colnames(pine.vst),1,1)),colSums)

heatmap.2(mean.dat[!is.na(rowSums(mean.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal)

#' #### Gene of interest (including Dark, sorted by tissue spec, 
#' mean of replicates)
heatmap.2(mean.dat[ord,][!is.na(rowSums(mean.dat[ord,])),],
          labRow=NA, dRowv = NULL, dendrogram = "column",
          trace="none",
          las=2,col=hpal)

#' #### Gene of interest (excluding Dark)
nd.sel <- substr(colnames(pine.vst),1,1)!="D"
s.dat <- s.dat[,nd.sel]

message(sprintf("There are %s genes of interest that are not expressed",
                sum(is.na(rowSums(s.dat)))))

heatmap.2(s.dat[!is.na(rowSums(s.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark, sorted by tissue spec)
exp.peak <- c("A","B")[apply(pESnD[rownames(s.dat),c("A","B")] == pESnD[rownames(s.dat),"maxn"],1,which)]
ord <- order(exp.peak,(1-pESnD[rownames(s.dat),"score"]))

heatmap.2(s.dat[ord,][!is.na(rowSums(s.dat[ord,])),],
          labRow=NA, dRowv = NULL, dendrogram = "column",
          trace="none",
          las=2,col=hpal,
          labCol = sub("\\.P\\.67N\\.","_",colnames(s.dat)))

#' #### Gene of interest (including Dark, mean of replicates)
mean.dat <- mean.dat[,colnames(mean.dat)!="D"]

heatmap.2(mean.dat[!is.na(rowSums(mean.dat)),],labRow=NA,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D")},
          trace="none",
          las=2,col=hpal)

#' #### Gene of interest (including Dark, sorted by tissue spec, 
#' mean of replicates)
heatmap.2(mean.dat[ord,][!is.na(rowSums(mean.dat[ord,])),],
          labRow=NA, dRowv = NULL, dendrogram = "column",
          trace="none",
          las=2,col=hpal)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
