#' ---
#' title: "P. taeda - create the STAR genome"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Aim
#' The goal is to create the STAR genome, i.e. scaffold contigs to make it
#' manageable in size. However, we want to avoid scaffolding read containing
#' regions, so we use an iteratively merged BAM file to avoid that. The BAM
#' file contains the STAR alignment of the digitally normalized catalog, ... DEFINE ME!
#' # Setup
#' ## Environment
#' Set the working dir
setwd("/mnt/picea/storage/reference/Pinus-taeda/v1.01/fasta")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/storage/reference/Pinus-taeda/v1.01/fasta")
#' ```
#' Load libraries
suppressPackageStartupMessages(library(Biostrings))

#' # Process
#' ## Read in the genomes
#' First the complete
ref <- readDNAStringSet("ptaeda.v1.01.scaffolds.fasta.gz")

#' Load libraries
suppressPackageStartupMessages(library(genomeIntervals))

#' Then the gene-containing only
#' This is for spruce:
#refGene <- readDNAStringSet("Pabies1.0-genome-gene-only.fa")

gff3 <- readGff3("../gff/Ptaeda1.01-Annotation-v3.0.gff3")

#' Sanity check
all(levels(seqnames(gff3)) %in% names(ref))

#' We select for the non-gene IDs
sel <- names(ref) %in% levels(seqnames(gff3))
refGene <- ref[sel]
ref <- ref[!sel]

#' And further select for size
sel <- width(ref) > 4999
rb5k <- ref[sel]
rs5k <- ref[!sel]

#' ## Scaffold
#' We build chromosomes  gapped by 10N using 
#' non gene scaffold smaller than 5k into 500Mb artefactual chromosomes
lim <- 514000000
pan <- DNAStringSet(rep(paste(rep("N",10),collapse=""),length(rs5k) * 2 -1))
pan[seq(1,length(pan),2)] <- rs5k

#' We also record the names
nams <- rep("spacer",length(rs5k) * 2 -1)
nams[seq(1,length(pan),2)] <- names(rs5k)

#' Before creating the chinks
max.chunks <- ceiling(sum(as.numeric(width(pan))) / lim)

#' We then define the position at which to split
starts <- match(0:(max.chunks-1),floor(cumsum(as.numeric(width(pan))) / lim))
starts <- cbind(starts,c(starts[-1]-1,length(pan)))
starts[,1] <- ifelse(starts[,1]%%2==0,starts[,1]+1,starts[,1])
starts[,2] <- ifelse(starts[,2]%%2==0,starts[,2]-1,starts[,2])

#' And finally construct the scaffold sequences. We also report the 
#' coordinates of the original scaffolds as a matrix
starts <- as.data.frame(starts)
starts$nam <- sprintf("Combined_%02d",1:max.chunks)

res <- apply(starts,1,function(ro,pan,nams){
  s <- as.integer(ro[1])
  e <- as.integer(ro[2])
  list(seq=DNAStringSet(unlist(pan[s:e])),
       nam=data.frame(
         scf=rep(ro[3],e-s+1),
         start=c(1,(cumsum(width(pan[s:e]))[-(e-s+1)]+1)),
         end=cumsum(width(pan[s:e])),
         nam=nams[s:e])
  )
},pan,nams)

#' Combine the sequences
seqs <- Reduce(append,lapply(res,"[[","seq"))
names(seqs) <- starts$nam

#' Simplify the names for STAR
names(refGene) <- sapply(strsplit(names(refGene)," "),"[",1)
names(rb5k) <- sapply(strsplit(names(rb5k)," "),"[",1)

#' Merge all sequences in a ref
ref <- c(refGene,rb5k,seqs)

#' And write it out
writeXStringSet(ref,file="ptaeda.v1.01-genome-collapsed-for-STAR.fa")


#' Then combine the mapping
mat <- do.call(rbind,lapply(res,"[[","nam"))

#' And write it out
write.table(mat,file="../annotation/ptaeda.v1.01-genome-collapsed-for-STAR-scaffold-coordinates.tsv",
            sep="\t",row.names=FALSE,quote=FALSE)

#' # Session Info
#' ```{r sessionInfo, echo=FALSE}
#' sessionInfo()
#' ```

