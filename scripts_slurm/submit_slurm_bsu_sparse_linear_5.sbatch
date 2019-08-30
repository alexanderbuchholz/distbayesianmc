#!/bin/bash
#
#SBATCH -J 5_HLA_FULL
#

#SBATCH -n 1              # on one node
#SBATCH -t 24-00:00          # Running time of 2 hours 24-0:00 (days-hours:minutes)


## 5 splits, 20 repetitions

#SBATCH --array=1-100%20

#SBATCH --mail-type=FAIL
#SBATCH -o /scratch/alexander/slurm_output/slurm-%A_%a.out # STDOUT

export TMPDIR=/scratch/alexander/distbayesianmc_sparselinear/tmp/

Rscript ../scripts_r/rscript_hla_sparse_ind.R $SLURM_ARRAY_TASK_ID "5" "split_loop"

