#' ---
#' title: "Sonali's Scots Pine light experiment Biological QA"
#' author: "Sonali Ranade and Nicolas Delhomme"
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
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Create a palette
#'pal = c("#00FFFF","#0000FF","#8A2BE2","#A52A2A","#DEB887","#5F9EA0","#7FFF00","#FF7F50","#DC143C","#00008B","#006400","#8B008B","#8B0000","#E9967A","#FFD700","#DAA520")
pal <- brewer.pal(8,"Dark2")
#'pal <- brewer.pal(12,"Set3")

#' Register the default plot margin
mar <- par("mar")

#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/spruce-pine-light/doc/pine.csv")

#' Read the HTSeq files in a matrix
res <- mclapply(dir("htseq-pine",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

names(res) <- sub(".*_P","P",sub("_sortmerna.*\\.txt","",dir("htseq-pine",pattern="*.txt")))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
samples <- samples[match(names(res),samples$SampleName),]
samples$SampleDescription <- factor(as.character(samples$SampleDescription))
samples$SampleName <- factor(as.character(samples$SampleName))

#' Raw Data QC analysis
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]
dir.create(file.path("analysis","HTSeq"),recursive = TRUE, showWarnings = FALSE)
write.csv(count.table,"analysis/HTSeq/raw-unormalised-data.csv")

#' Extract the HTSeq stat lines
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' Convert them into percentages
pander(apply(count.stats,2,function(co){round(co*100/sum(co))}))

#' Plot the stats
#' 
#' There are no outliers
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6),cex.names=.6)
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

#' The average percentage of aligned reads is 57%
# round(mean(unlist(count.stats["aligned",]/colSums(count.stats))),digits=2)*100
boxplot(unlist(count.stats["aligned",]/colSums(count.stats)),
        main="aligned reads",ylab="percent aligned",ylim=c(0,1))

#' Check how many genes are never expressed
sel <- rowSums(count.table) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count.table),digits=1),
        sum(sel),
        nrow(count.table))

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative coverage is as expected, around 100X
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by sample type
#' 
plot.multidensity(log10(count.table),
                  col=c(1,2,pal)[as.integer(samples$SampleDescription)],
                  legend.x="topright",
                  legend=levels(samples$SampleDescription),
                  legend.col=c(1,2,pal),
                  legend.lwd=2,legend.cex = 0.7,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' # Data normalisation 
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' Create the dds object
conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' Perform the VST
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)
vst <- vst - min(vst)
write.csv(vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_data.csv")

#' Validate the VST 
#' 
#' Visualize the corrected mean - sd relationship. It is fairly linear,
#' meaning we can assume homoscedasticity.
#' The slight initial trend / bump is due to genes having few counts in
#' a few subset of the samples and hence having a higher variability. This is
#' expected.
meanSdPlot(vst[rowSums(count.table)>0,])

#' # QC on the normalised data
#' 
#' ## PCA
#' 
#' First perform a Principal Component Analysis (PCA) of the data
#'to do a quick quality assessment; i.e. replicate should cluster
#' and the first 2-3 dimensions shouldbe explainable by biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Plot the PCA 3 first dimensions
#' There are too many points to be informative. There seem to be 3 groups but
#' possibly 2 sets of outliers
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=rep(c(1,2,3,pal),2)[as.integer(samples$SampleDescription)],
              pch=19)
par(mar=mar)

#' Then the first two dimensions
#' The first dimension separates the light conditions (Dark, A & B). 
#' The second dimension shows variability within the two ecotypes
#' The Dark samples are more spread out than the light ones. The light ones
#' do not show a large ecotype variance. A number of P56 samples under the
#' A light condition display appear to be outlier. In a more limited effect
#' a number of P67 samples also appear to be outliers in the Dark conditions.
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=c(pal[1:3])[as.integer(factor(substr(samples$SampleDescription,1,1)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("top",col=pal[1:3],
       legend=levels(factor(substr(samples$SampleDescription,1,1))),
       pch=19,bty="n")

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=c(pal[1:2])[as.integer(factor(substr(samples$SampleDescription,3,5)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("top",col=pal[1:2],
       legend=levels(factor(substr(samples$SampleDescription,3,5))),
       pch=19,bty="n")


#' And the 2nd and 3rd dims
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=c(pal[1:3])[as.integer(factor(substr(samples$SampleDescription,1,1)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("top",col=pal[1:3],
       legend=levels(factor(substr(samples$SampleDescription,1,1))),
       pch=19,bty="n")

plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=c(pal[1:2])[as.integer(factor(substr(samples$SampleDescription,3,5)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("top",col=pal[1:2],
       legend=levels(factor(substr(samples$SampleDescription,3,5))),
       pch=19,bty="n")

#' Identify these samples
pander(as.character(unique(samples$SampleDescription[which(pc$x[,2] < -40)])))

#' ## Heatmap
#' The 1000 most variable genes are selected and plotted as a heatmap
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6,
          labCol =  samples$SampleDescription)

#' The distinction between Dark and Light treated (A & B) is very clear 
#' in these 1000 first genes. The technical replicates are very similar. The
#' outliers A-P56-2 and D-P67-1 are closer to the Dark part of the cluster. The
#' ecotype difference is not obvious though apart from the Dark samples - if 
#' we consider D-P67-1 an outlier.
#' ```{r empty, eval=FALSE, echo=FALSE}
#' ```
#' # Conclusion
#' The raw quality of the data appears very good. The data normalisation
#' gives satisfying results (as far as the VST fit is concerned). The PCA and 
#' the heatmap identifies two samples as possible outliers. Moreover, it is clear
#' from the analysis that the technical replicates can be combined together,
#' which is what we do next

count.table <- do.call(cbind,lapply(split.data.frame(t(count.table),samples$SampleName),colSums))
colnames(count.table) <- samples$SampleDescription[match(colnames(count.table),samples$SampleName)]
count.table <- count.table[,order(colnames(count.table))]
write.csv(count.table,"analysis/HTSeq/raw-unormalised_tech-rep-combined_data.csv")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
