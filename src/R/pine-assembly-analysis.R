#' ---
#' title: "Pine de-novo transcriptome assembly analysis"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/facility/gmap")

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/facility/gmap")
#' ```

#' Load libraries
# suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(genomeIntervals))
suppressPackageStartupMessages(library(GenomicRanges))
#suppressPackageStartupMessages(library(LSD))
#suppressPackageStartupMessages(library(pander))
#suppressPackageStartupMessages(library(RColorBrewer))
# suppressPackageStartupMessages(library(VennDiagram))

#' Helpers
# source("~/Git/UPSCb/src/R/percentile.R")
source("~/Git/UPSCb/src/R/gff3Utilities.R")
#source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Graphic parameters
#mar <- par("mar")
#pal <- brewer.pal(8,"Dark2")

#' # Data
#' ## Annot
gff <- readGff3("/mnt/picea/storage/reference/Pinus-taeda/v1.01/gff/Ptaeda1.01-Annotation-v3.0_synthetic-transcripts.gff3",
                       quiet=TRUE)

#' Select the genes only
gff <- gff[gff$type == "gene",]
grng <- as(gff,"GRanges")

#' ## GMAP results
GMAP.files <- dir(pattern="\\.gff3\\.[m,t,u].*",full.names=TRUE)
gff_GMAPs <- lapply(GMAP.files,
                    readGff3,
                    quiet=TRUE)

names(gff_GMAPs) <- sub(".*\\.","",GMAP.files)
gff_GMAPs <- lapply(gff_GMAPs,function(g){g[g$type=="gene"]})
grng_GMAPS <- lapply(gff_GMAPs,as,"GRanges")

#' Number of alignments
elementNROWS(grng_GMAPS)

#' Overlaps
ovls <- lapply(grng_GMAPS,findOverlaps,grng)

#' number of Overlaps
elementNROWS(ovls)

length(unique(queryHits(ovls$uniq)))
length(unique(subjectHits(ovls$uniq)))

length(unique(queryHits(ovls$transloc)))
length(unique(subjectHits(ovls$transloc)))

length(unique(queryHits(ovls$mult)))
length(unique(subjectHits(ovls$mult)))

length(unique(unlist(lapply(lapply(ovls,subjectHits),unique))))
length(unique(unlist(lapply(lapply(ovls,queryHits),unique))))

res <- lapply(1:3,function(i,qlist){
    reduce(grng_GMAPS[[i]][qlist[[i]],])
},lapply(ovls,queryHits))

reduce(Reduce("c",res))

20099 / 20361



table(1:length(grng_GMAPS[["uniq"]]) %in% unique(queryHits(ovls[["uniq"]])))

boxplot(width(grng_GMAPS[["uniq"]][! 1:length(grng_GMAPS[["uniq"]]) %in% unique(queryHits(ovls[["uniq"]])),]),
        width(grng_GMAPS[["uniq"]][1:length(grng_GMAPS[["uniq"]]) %in% unique(queryHits(ovls[["uniq"]])),]),
        names=c("not","aligned"),log="y",)

t.test(width(grng_GMAPS[["uniq"]][! 1:length(grng_GMAPS[["uniq"]]) %in% unique(queryHits(ovls[["uniq"]])),]),
       width(grng_GMAPS[["uniq"]][1:length(grng_GMAPS[["uniq"]]) %in% unique(queryHits(ovls[["uniq"]])),]))

median(width(grng_GMAPS[["uniq"]][! 1:length(grng_GMAPS[["uniq"]]) %in% unique(queryHits(ovls[["uniq"]])),]))
median(width(grng_GMAPS[["uniq"]][ 1:length(grng_GMAPS[["uniq"]]) %in% unique(queryHits(ovls[["uniq"]])),]))

reduce(grng_GMAPS[["uniq"]][! 1:length(grng_GMAPS[["uniq"]]) %in% unique(queryHits(ovls[["uniq"]])),])

table(1:length(grng_GMAPS[["transloc"]]) %in% unique(queryHits(ovls[["transloc"]])))
table(1:length(grng_GMAPS[["mult"]]) %in% unique(queryHits(ovls[["mult"]])))
