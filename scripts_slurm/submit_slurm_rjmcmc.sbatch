#!/bin/bash
#
#SBATCH -J rjmcm
#

#SBATCH -n 1              # on one node
#SBATCH -t 12-0:00         # Running time of 2 hours

#SBATCH --array=1-20


Rscript ../scripts_r/rscript_rjmcmc.R $SLURM_ARRAY_TASK_ID

