#!/bin/bash
#
#SBATCH -J flightcomplex
#

#SBATCH -n 1              # on one node
#SBATCH -t 12-0:00         # Running time of 2 hours

#SBATCH --array=1-40


Rscript rscript_modelsplitting.R "30" $SLURM_ARRAY_TASK_ID "params_logistic.R" "flights_complex2"
Rscript rscript_modelsplitting.R "50" $SLURM_ARRAY_TASK_ID "params_logistic.R" "flights_complex2"
Rscript rscript_modelsplitting.R "100" $SLURM_ARRAY_TASK_ID "params_logistic.R" "flights_complex2"
Rscript rscript_modelsplitting.R "200" $SLURM_ARRAY_TASK_ID "params_logistic.R" "flights_complex2"

Rscript rscript_modelsplitting.R "30" $SLURM_ARRAY_TASK_ID "params_logistic.R" "flights_complex1"
Rscript rscript_modelsplitting.R "50" $SLURM_ARRAY_TASK_ID "params_logistic.R" "flights_complex1"
Rscript rscript_modelsplitting.R "100" $SLURM_ARRAY_TASK_ID "params_logistic.R" "flights_complex1"
Rscript rscript_modelsplitting.R "200" $SLURM_ARRAY_TASK_ID "params_logistic.R" "flights_complex1"
