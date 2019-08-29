#!/bin/bash
#
#SBATCH -J 1000_HIGGS_FULL
#

#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH -p skylake-himem
#SBATCH -n 1              # on one node
#SBATCH -t 24:00:00          # Running time of 2 hours 24-0:00 (days-hours:minutes)

## 1000 splits, 20 repetitions

#SBATCH --array=1-1000

#SBATCH --mail-type=FAIL
#SBATCH -o /mrc-bsu/scratch/ab2603/slurm_output/slurm-%A_%a.out # STDOUT


export TMPDIR=/mrc-bsu/scratch/ab2603/distbayesianmc_higgs/tmp/

module load r-3.6.0-gcc-5.4.0-bzuuksv  
Rscript ../scripts_r/rscript_HIGGS_logit_ind.R $SLURM_ARRAY_TASK_ID "1000" "full_loop"

