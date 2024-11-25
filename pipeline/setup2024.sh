#!/bin/bash -l
set -ex
proj=/mnt/picea/projects/spruce/facility/

# Raw data
mkdir data
ln -s $proj/rgarcia-pine/raw data/pine
ln -s $proj/rgarcia-spruce/raw data/spruce

# Analysis dir
mkdir analysis
ln -s $proj/rgarcia-pine/analysis analysis/pine
ln -s $proj/rgarcia-spruce/analysis analysis/spruce
