#!/bin/bash
#
#SBATCH -J 100_HIGGS_FULL
#

#SBATCH -n 1              # on one node
#SBATCH -t 24-00:00          # Running time of 2 hours 24-0:00 (days-hours:minutes)

## 100 splits, 20 repetitions

#SBATCH --array=1-2000%40 

#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/alexander/slurm_output/slurm-%A_%a.out # STDOUT

export TMPDIR=/scratch/alexander/distbayesianmc_higgs/tmp/

Rscript ../scripts_r/rscript_HIGGS_logit_ind.R $SLURM_ARRAY_TASK_ID "100" "split_loop"

