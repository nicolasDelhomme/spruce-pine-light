#!/bin/bash -l

set -ex

mail="sonali.ranade@slu.se"
proj=b2015329
in=/proj/$proj/nobackup/pine-two

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

sbatch -A $proj -e $in/reports.err -o $in/reports.out --mail-user $mail $UPSCb/pipeline/runRnaSeqPipelineStats.sh $in
