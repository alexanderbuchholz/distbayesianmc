#!/bin/bash
#
#SBATCH -J 100_HIGGS_FULL
#

#SBATCH -n 1              # on one node
#SBATCH -t 24-00:00          # Running time of 2 hours 24-0:00 (days-hours:minutes)

#SBATCH --array=1-2000%40

#SBATCH --mail-type=FAIL

export TMPDIR=/scratch/alexander/$SLURM_JOBID

Rscript ../scripts_r/rscript_HIGGS_logit_ind.R $SLURM_ARRAY_TASK_ID "100"

