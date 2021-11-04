#' ---
#' title: "Convert the P. tremula gff3 to a format PASA understands"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/storage/reference/Pinus-taeda/v1.01/gff")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/storage/reference/Pinus-taeda/v1.01/gff")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(genomeIntervals))
suppressPackageStartupMessages(library(GenomicRanges))

#' # Process
#' Read the gff file
gff3 <- readGff3("Ptaeda1.01-Annotation-v3.0.gff3")

#' Extract the exons
exons <- gff3[gff3$type=="exon",]
exonID <- getGffAttribute(exons,"ID")
mRnaID <- strsplit(getGffAttribute(exons,"Parent"),",")

#' Check for duplication
table( elementLengths(mRnaID) )
grngs <- as(exons,"GRanges")
length(grngs)
ovl <- findOverlaps(grngs,type="equal",ignoreSelf=TRUE,ignoreRedundant=TRUE)

#' Exons have been duplicated already, no need to do it