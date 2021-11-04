#!/bin/bash -l

## error and verbose
set -ex

## input
inxDir=/mnt/picea/storage/reference/Pinus-taeda/v1.01/indices/STAR
fasta=/mnt/picea/storage/reference/Pinus-taeda/v1.01/fasta/ptaeda.v1.01-genome-collapsed-for-STAR.fa
gtf=/mnt/picea/storage/reference/Pinus-taeda/v1.01/gtf/Ptaeda1.01-Annotation-v3.0.gtf

## prepare
mkdir -p $inxDir/Ptaeda1.01-genome
sbatch -w watson -e $inxDir/Ptaeda1.01-genome.err -o $inxDir/Ptaeda1.01-genome.out --mem=500G -n 64 --mail-user="nicolas.delhomme@umu.se" $UPSCb/pipeline/runSTARGenomeGenerate.sh -m 500000000000 -p 64 -f gtf -b 17 $inxDir/Ptaeda1.01-genome $fasta $gtf
