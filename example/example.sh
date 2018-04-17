#!/bin/bash
#SBATCH --job-name=example                                           # job name
#SBATCH --partition=super                                       # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                                    # number of nodes requested by user
#SBATCH --time=10-00:00:00                                            # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=./sbatch_output_%j                                 # redirect both standard output and erro output to the same file
#SBATCH --error=./sbatch_error_%j
source /home2/twang6/.bash_profile

###########JOBSTART########################

echo "e"
