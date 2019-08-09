#!/bin/bash
#
#SBATCH -p skylake # Partition normal
#SBATCH -J rjmcmc
#

#SBATCH -n 1              # on one node
#SBATCH -t 6-0:00         # Running time of 2 hours

#SBATCH --array=1-20

. /etc/profile.d/modules.sh                # This line enables the module command
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
module load r-3.6.0-gcc-5.4.0-bzuuksv 
module load jdk-8u131-b11-gcc-5.4.0-z7brlnq


Rscript ../scripts_r/rscript_rjmcmc.R $SLURM_ARRAY_TASK_ID



