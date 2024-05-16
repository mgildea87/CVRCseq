#!/bin/bash

# Setting environment variables
export TMPROOT=/gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq # environment path
export PATH=$TMPROOT/bin:$PATH
export PYTHONPATH=$TMPROOT/lib/python3.10/site-packages:$PYTHONPATH # replace * with the python version being used
#Custom packages should be set in PYTHONPATH as well
export LD_LIBRARY_PATH=$TMPROOT/lib:$LD_LIBRARY_PATH

# Activating environment
module purge
module load slurm
source /gpfs/share/apps/anaconda3/gpu/new/etc/profile.d/conda.sh #should be replaced with conda.sh of anaconda module environment was created with
conda activate /gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq # may need full path
