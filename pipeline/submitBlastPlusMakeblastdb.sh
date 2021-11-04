#!/bin/bash -l

## error and verbose
set -ex

## input
inxDir=/mnt/picea/storage/reference/Pinus-taeda/v1.01/indices/BLAST+
fasta=/mnt/picea/storage/reference/Pinus-taeda/v1.01/fasta/ptaeda.v1.01.scaffolds.trimmed.all.maker.proteins.fa
title=Ptaeda1.01-peptide

## prepare
mkdir -p $inxDir
sbatch -e $inxDir/Ptaeda1.01-peptide.err -o $inxDir/Ptaeda1.01-peptide.out \
--mail-user="nicolas.delhomme@umu.se" $UPSCb/pipeline/runBlastPlusMakeblastdb.sh \
-p prot -t $title $fasta $inxDir
