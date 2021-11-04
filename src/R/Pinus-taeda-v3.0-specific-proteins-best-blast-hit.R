#' ---
#' title: "Sari's Scots Pine v3.0 specific proteins best blast hits"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/storage/reference/Pinus-taeda/v1.01/annotation/BLAST+")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/storage/reference/Pinus-taeda/v1.01"/annotation/BLAST+)
#' ```

#' Helper functions
source("~/Git/UPSCb/src/R/blastUtilities.R")

#' VARS
format <- c("query.id",
            "subject.id",
            "percent.identity",
            "alignment.length",
            "mismatches",
            "gap.opening",
            "query.start",
            "query.end",
            "subject.start",
            "subject.end",
            "e.value",
            "bit.score",
            "query.length",
            "subject.length")

#' # Process
#' Read the blt
blt <- readBlast("Pabies1.0-all.phase.gff3.AA.fa_Ptaeda1.01-Annotation-v3.0-specific-proteins.blt.gz",
                 format,bestHit = TRUE)

#' Extract the best blast hits
write.table(blt$df,sep="\t",quote = FALSE,
            file="Pabies1.0_Ptaeda1.01-Annotation-v3.0-specific-proteins_best-BLAST-hit.tsv")

#' Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```


