#!/bin/bash -l

set -ex

proj=b2015329
mail=sonali.ranade@slu.se
in=/proj/$proj/nobackup/pine-two/raw
out=/proj/$proj/nobackup/pine-two
genome=/proj/b2011227/indices/STAR/Ptaeda1.01-genome
gff3=/proj/b2011227/reference/gff3/Ptaeda1.01-Annotation-v3.0.gff3
hgff3=/proj/b2011227/reference/gff3/Ptaeda1.01-Annotation-v3.0_synthetic-transcripts.gff3
start=1
end=8
mem=512

module load bioinfo-tools samtools star/2.4.0f1 htseq

export SORTMERNADIR=/home/delhomme/sortmerna

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fastq.gz"`; 
do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
  bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -g $genome \
  -G $gff3 -H $hgff3 -t -m $mem $proj $mail ${line}_1.fastq.gz ${line}_2.fastq.gz $out
done
