#!/bin/bash
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH -A kaw

cd /mnt/picea/home/fai/git/spruce-pine-light
nextflow run nf-core/rnaseq -r 3.14.0 -profile upscb,singularity -work-dir analysis/pine/work -params-file nextflow/nf-params-pine.json -c nextflow/upscb.config
