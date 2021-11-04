#!/bin/bash -l

## error and verbose
set -ex

## check for the UPSCb env. var.
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
    usage
fi

## default args
#in=/mnt/picea/storage/reference/Pinus-taeda/v1.01/fasta/ptaeda.v1.01.scaffolds.trimmed.all.maker.proteins.fa
in=/mnt/picea/storage/reference/Pinus-taeda/v1.01/fasta/Ptaeda1.01-Annotation-v3.0-specific-proteins.fa
out=/mnt/picea/storage/reference/Pinus-taeda/v1.01/annotation/BLAST+
#inxDir=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/BLAST+
inxDir=/mnt/picea/storage/reference/UniRef90/201701/indices/blast+2.6.0
#inx=Pabies1.0-all.phase.gff3.AA.fa
inx=uniref90
CPU=60

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch -e $out/Ptaeda1.01-peptide.err -o $out/Ptaeda1.01-peptide.out \
-c $CPU $UPSCb/pipeline/runBlastPlus.sh -p $CPU blastp $in $inxDir/$inx $out
