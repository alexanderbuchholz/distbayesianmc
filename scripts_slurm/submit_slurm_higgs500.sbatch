#!/bin/bash
#
#SBATCH -J 500_HIGGS_FULL
#

#SBATCH -n 1              # on one node
#SBATCH -t 24-00:00          # Running time of 2 hours 24-0:00 (days-hours:minutes)

## 500 splits, 20 repetitions

#SBATCH --array=1-10000%40

#SBATCH --mail-type=FAIL

export TMPDIR=/scratch/alexander/distbayesianmc_higgs/tmp/

Rscript ../scripts_r/rscript_HIGGS_logit_ind.R $SLURM_ARRAY_TASK_ID "500"

