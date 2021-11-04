#!/bin/bash -l
set -ex
proj=b2015329

# Pine
## create the dirs 
mkdir -p /proj/$proj/nobackup/pine-light/raw
cd /proj/$proj/nobackup/pine-light/raw

ln -s /proj/$proj/sequence-data/R.GarciaGil_15_01/P2254_1[0-2][0-9]/*/*.fastq.gz .
ln -s /proj/$proj/sequence-data/R.GarciaGil_15_01/P2254_13[0-6]/*/*.fastq.gz .


# Spruce
## create the dirs
mkdir -p /proj/$proj/nobackup/spruce-light/raw
cd /proj/$proj/nobackup/spruce-light/raw

ln -s /proj/$proj/sequence-data/R.GarciaGil_15_01/P2254_13[8,9]/*/*.fastq.gz .
ln -s /proj/$proj/sequence-data/R.GarciaGil_15_01/P2254_1[^0-3][0-9]/*/*.fastq.gz .
