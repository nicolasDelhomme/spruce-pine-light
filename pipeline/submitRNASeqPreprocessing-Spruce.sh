#!/bin/bash -l

set -ex

proj=b2015329
mail=sonali.ranade@slu.se
## in=/proj/$proj/nobackup/spruce-light/raw
## out=/proj/$proj/nobackup/spruce-light
in=/proj/$proj/nobackup/spruce-two/raw
out=/proj/$proj/nobackup/spruce-two
genome=/proj/b2011227/indices/STAR/Pabies01-genome
gff3=/proj/b2011227/reference/gff3/Eugene.gff3
start=3
end=8
mem=fat

module load bioinfo-tools samtools star/2.4.0f1 htseq

export SORTMERNADIR=/home/delhomme/sortmerna

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fastq.gz"`; 
do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
  bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -g $genome \
  -G $gff3 -H $gff3 -t -m $mem $proj $mail ${line}_1.fastq.gz ${line}_2.fastq.gz $out
done
