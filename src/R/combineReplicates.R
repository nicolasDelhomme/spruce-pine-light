#' ---
#' title: "Sari's Norway Spruce light experiment - technical replicates merge "
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

#' Libraries
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))

#' A palette
pal <- brewer.pal(3,"Dark2")

#' # Process
#' Read the sample information
spruce.samples <- read.csv("~/Git/UPSCb/projects/spruce-pine-light/doc/spruce.csv")
pine.samples <- read.csv("~/Git/UPSCb/projects/spruce-pine-light/doc/pine.csv")

#' Read the spruce HTSeq files in a list
spruce.res <- mclapply(dir("htseq-spruce",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

#' And the same for the pine
pine.res <- mclapply(dir("htseq-pine",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

#' Rename the samples
names(spruce.res) <- spruce.samples[match(sub(".*_P","P",
                        sub("_sortmerna.*\\.txt","",
                            dir("htseq-spruce",pattern="*.txt"))),
                    spruce.samples$SampleName),"SampleDescription"]

names(pine.res) <- pine.samples[match(sub(".*_P","P",
                                          sub("_sortmerna.*\\.txt","",
                                              dir("htseq-pine",pattern="*.txt"))),
                                      pine.samples$SampleName),"SampleDescription"]

#' Merge the samples (2 for spruce, 4 for pine)
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual",
             "__not_aligned","__alignment_not_unique")

#' Spruce
sel <- match(addInfo,spruce.res[[1]][,1])
spruce.mat <- sapply(lapply(unique(names(spruce.res)),grep,names(spruce.res)),
       function(pos,res,sel){
         if(length(pos)>1){
           cts <- do.call("+",lapply(res[pos[1:2]],"[",3))[-sel,1]
         } else {
           cts <- res[[pos]][-sel,3]
         }
         names(cts) <- res[[pos[1]]][-sel,1]
         return(cts)
       },spruce.res,sel)
colnames(spruce.mat) <- unique(names(spruce.res))

#' Pine
sel <- match(addInfo,pine.res[[1]][,1])
pine.mat <- sapply(lapply(unique(names(pine.res)),grep,names(pine.res)),
                     function(pos,res,sel){
                       if(length(pos)>1){
                         cts <- rowSums(do.call(cbind,sapply(res[pos[1:4]],"[",2)))[-sel]
                       } else {
                         cts <- res[[pos]][-sel,2]
                       }
                       names(cts) <- res[[pos[1]]][-sel,1]
                       return(cts)
                     },pine.res,sel)
colnames(pine.mat) <- unique(names(pine.res))

#' # Validate
#' This is just a validation on the raw data - just to make sure we have not 
#' wrongly combined samples

#' ## Sequencing depth
barplot(colSums(spruce.mat),main="sequencing depth",las=2)
barplot(colSums(pine.mat),main="sequencing depth",las=2)

#' ## PCA
#' Spruce
pc <- prcomp(t(spruce.mat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(spruce.mat),1,1)))],
     pch=c(17,19)[as.integer(factor(substr(colnames(spruce.mat),4,5)))],
     main="Principal Component Analysis",sub="variance stabilized counts")

#' Pine
pc <- prcomp(t(pine.mat))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(substr(colnames(pine.mat),1,1)))],
     pch=c(17,19)[as.integer(factor(substr(colnames(pine.mat),4,5)))],
     main="Principal Component Analysis",sub="variance stabilized counts")

#' ## Heatmap
#' Spruce
sel <- order(apply(spruce.mat,1,sd),decreasing=TRUE)[1:1000]
heatmap.2(log10(spruce.mat[sel,]+1),labRow = NA,trace = "none",cexCol = 0.6)

#' Pine
sel <- order(apply(pine.mat,1,sd),decreasing=TRUE)[1:1000]
heatmap.2(log10(pine.mat[sel,]+1),labRow = NA,trace = "none",cexCol = 0.6)

#' # Export
#' Spruce
write.csv(spruce.mat,"analysis/HTSeq/spruce_raw-unormalised_technical-replicates-merged_data.csv")

#' Pine
write.csv(pine.mat,"analysis/HTSeq/pine_raw-unormalised_technical-replicates-merged_data.csv")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

