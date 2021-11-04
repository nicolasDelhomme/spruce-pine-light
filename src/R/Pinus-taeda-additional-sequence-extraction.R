#' ---
#' title: "Sari's Scots Pine transcript annotation"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/storage/reference/Pinus-taeda/v1.01")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/storage/reference/Pinus-taeda/v1.01")
#' ```

#' Libraries
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(genomeIntervals))
suppressPackageStartupMessages(library(VennDiagram))

#' Helper functions
source("~/Git/UPSCb/src/R/sequenceUtilities.R")

#' #  Process
#' ## Extraction
#' ### Read the complete annotation
gff3 <- readGff3("gff/Ptaeda1.01-Annotation-v3.0.gff3")

#' ### Read the genome sequence
genome <- readDNAStringSet("fasta/ptaeda.v1.01.scaffolds.fasta.gz")

#' ### Extract the transcripts
transcritps <- extractMrnaFromGenome(gff3,genome)

#' ### Extract the proteins
proteins <- extractProteinFromGenome(gff3,genome)
pstate <- getProteinState(proteins)
barplot(table(pstate))
barplot(table(pstate),log="y")

#' ### Validation
#' #### Comparison with the existing sequences
known.proteins <- readAAStringSet("fasta/ptaeda.v1.01.scaffolds.trimmed.all.maker.proteins.fa")
names(known.proteins) <- sub(" .*","",names(known.proteins))
pkstate <- getProteinState(known.proteins)
barplot(table(pkstate))

# Overlap of the names
grid.draw(venn.diagram(list(
  "v3"=names(proteins),
  "v1.01"=names(known.proteins)),
    filename=NULL))

#' Find the common sequences (the majority)
nams <- intersect(names(proteins),names(known.proteins))

#' Drop the stop codons
no.stop.proteins <- AAStringSet(sub("\\*","",proteins[nams]))

#' Compare the sizes
table(width(no.stop.proteins) - width(known.proteins[nams]))

sel <- abs(width(no.stop.proteins) - width(known.proteins[nams])) == 36
no.stop.proteins[sel]
known.proteins[names(no.stop.proteins[sel])]

#' Two IDs have been inverted...

#' ## Export
#' ### All sequences
writeXStringSet(proteins,"fasta/Ptaeda1.01-Annotation-v3.0-proteins.fa")
writeXStringSet(transcritps,"fasta/Ptaeda1.01-Annotation-v3.0-transcripts.fa")

#' ### Extracting the novel sequences
writeXStringSet(proteins[setdiff(names(proteins),names(known.proteins))],
                "fasta/Ptaeda1.01-Annotation-v3.0-specific-proteins.fa")

writeXStringSet(known.proteins[setdiff(names(known.proteins),names(proteins))],
                "fasta/Ptaeda1.01-Annotation-v1.0-specific-proteins.fa")

write(sprintf("The IDs %s and %s are inverted between v1.0 and v3.0",
      names(no.stop.proteins[sel])[1],names(no.stop.proteins[sel])[2]),file = "WARNING")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
