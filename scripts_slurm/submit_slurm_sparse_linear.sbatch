#!/bin/bash
#
#SBATCH -J sparse_linear
#

#SBATCH -n 1              # on one node
#SBATCH -t 12-0:00         # Running time of 2 hours

#SBATCH --array=1-40


Rscript ../scripts_r/rscript_modelsplitting.R "5" $SLURM_ARRAY_TASK_ID "params_sparse" "sparse_reg_1"
Rscript ../scripts_r/rscript_modelsplitting.R "10" $SLURM_ARRAY_TASK_ID "params_sparse" "sparse_reg_1"
Rscript ../scripts_r/rscript_modelsplitting.R "20" $SLURM_ARRAY_TASK_ID "params_sparse" "sparse_reg_1"
Rscript ../scripts_r/rscript_modelsplitting.R "50" $SLURM_ARRAY_TASK_ID "params_sparse" "sparse_reg_1"

Rscript ../scripts_r/rscript_modelsplitting.R "5" $SLURM_ARRAY_TASK_ID "params_sparse" "sparse_reg_2"
Rscript ../scripts_r/rscript_modelsplitting.R "10" $SLURM_ARRAY_TASK_ID "params_sparse" "sparse_reg_2"
Rscript ../scripts_r/rscript_modelsplitting.R "20" $SLURM_ARRAY_TASK_ID "params_sparse" "sparse_reg_2"
Rscript ../scripts_r/rscript_modelsplitting.R "50" $SLURM_ARRAY_TASK_ID "params_sparse" "sparse_reg_2"
