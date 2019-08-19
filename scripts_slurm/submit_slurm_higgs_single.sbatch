#!/bin/bash
#
#SBATCH -J 100_HIGGS_FULL
#

#SBATCH -n 1              # on one node
#SBATCH -t 24-00:00          # Running time of 2 hours 24-0:00 (days-hours:minutes)

##SBATCH --array=1-100%40

#SBATCH --mail-type=FAIL

Rscript ../scripts_r/rscript_HIGGS_logit_ind.R "1" "1000"

