#' ---
#' title: "P. taeda - create the synthetic transcripts"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Aim
#' The goal is to create synthetic transcripts from the Pinus taeda v1.01 gff3 
#' file
#' # Setup
#' ## Environment
#' Set the working dir
setwd("/mnt/picea/storage/reference/Pinus-taeda/v1.01/gff")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/storage/reference/Pinus-taeda/v1.01/gff")
#' ```
#' Load libraries
suppressPackageStartupMessages(library(easyRNASeq))

#' # Process
gff3.filename <- "Ptaeda1.01-Annotation-v3.0.gff3"
gI <- createSyntheticTranscripts(gff3.filename,output="Genome_intervals")

#' # Export
writeGff3(gI,"Ptaeda1.01-Annotation-v3.0_synthetic-transcripts.gff3")

#' # SessionInfo
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
