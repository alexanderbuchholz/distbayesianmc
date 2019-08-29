#!/bin/bash
#
#SBATCH -J test_1000_HLA_FULL
#

#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH -p skylake-himem

#SBATCH -n 1              # on one node
#SBATCH -t 2:00:00          # Running time of 2 hours 24-0:00 (days-hours:minutes)

##SBATCH --array=1-100%40

#SBATCH --mail-type=FAIL
#SBATCH -o /mrc-bsu/scratch/ab2603/slurm_output/slurm-%A_%a.out # STDOUT

export TMPDIR=/scratch/alexander/distbayesianmc_sparselinear/tmp/
module load r-3.6.0-gcc-5.4.0-bzuuksv
Rscript ../scripts_r/rscript_hla_sparse_ind.R "1" "1000" "split_loop"

