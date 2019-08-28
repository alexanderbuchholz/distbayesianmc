#!/bin/bash
#
#SBATCH -J rjmcmc
#

#SBATCH -n 1              # on one node
#SBATCH -t 12-0:00         # Running time of 2 hours

#SBATCH --array=1-20

#SBATCH --mail-type=FAIL

#SBATCH -o /scratch/alexander/slurm_output/slurm-%A_%a.out # STDOUT
export TMPDIR=/scratch/alexander/distbayesianmc_rjmcmc/tmp/

Rscript ../scripts_r/rscript_rjmcmc.R $SLURM_ARRAY_TASK_ID

