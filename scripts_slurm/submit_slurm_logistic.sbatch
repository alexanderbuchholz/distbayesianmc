#!/bin/bash
#
#SBATCH -J logistic
#

#SBATCH -n 1              # on one node
#SBATCH -t 24-0:00         # Running time of 2 hours

#SBATCH --array=1-20

#SBATCH --mail-type=FAIL

Rscript ../scripts_r/rscript_logit_ind.R $SLURM_ARRAY_TASK_ID

