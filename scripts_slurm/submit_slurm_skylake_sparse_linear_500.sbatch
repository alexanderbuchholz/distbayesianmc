#!/bin/bash
#
#SBATCH -J 500_HLA_FULL
#

#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH -p skylake-himem
#SBATCH -n 1              # on one node
#SBATCH -t 6:00          # Running time of 2 hours 24-0:00 (days-hours:minutes)

## 500 splits, 20 repetitions

#SBATCH --array=1-10000%500

#SBATCH --mail-type=FAIL
#SBATCH -o /mrc-bsu/scratch/ab2603/slurm_output/output.%a.out # STDOUT

export TMPDIR=/mrc-bsu/scratch/ab2603/distbayesianmc_sparselinear/tmp/
  
  
module load r-3.6.0-gcc-5.4.0-bzuuksv
Rscript ../scripts_r/rscript_hla_sparse_ind.R $SLURM_ARRAY_TASK_ID "500"

