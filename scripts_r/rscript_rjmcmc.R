#!/usr/bin/env Rscript
library(distbayesianmc)
#task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#task_id <- as.numeric(task_id_string)
args = commandArgs(trailingOnly=TRUE)
print(args)

sim_id <- as.numeric(args[1])

list_params_model_onesplit <- list(scale =  1,
                                   ssplits =  1,
                                   typesplit =  "random",
                                   dataset =  "higgs1_small", 
                                   mrep = 24,
                                   n.mil = 5, 
                                   n.cores = 1)

list_params_model_multisplits1 <- list(scale =  1,
                                       ssplits =  2,
                                       typesplit =  "random",
                                       dataset =  "higgs1_small",#"",
                                       mrep = 24, 
                                       n.mil = 5,
                                       n.cores = 1)

list_params_model_multisplits2 <- list(scale =  1,
                                       ssplits =  3,
                                       typesplit =  "random",
                                       dataset =  "higgs1_small",#"",
                                       mrep = 24,
                                       n.mil = 5,
                                       n.cores = 1)

list_params_model_multisplits3 <- list(scale =  1,
                                       ssplits =  5,
                                       typesplit =  "random",
                                       dataset =  "higgs1_small",#"",
                                       mrep = 24,
                                       n.mil = 5,
                                       n.cores = 1)

list_params_model_multisplits4 <- list(scale =  1,
                                       ssplits =  10,
                                       typesplit =  "random",
                                       dataset =  "higgs1_small",#"",
                                       mrep = 24,
                                       n.mil = 5,
                                       n.cores = 1)


list_params_model <- list(list_params_model_onesplit, 
                          list_params_model_multisplits1, 
                          list_params_model_multisplits2, 
                          list_params_model_multisplits3, 
                          list_params_model_multisplits4)

setwd("/scratch/alexander/distbayesianmc_rjmcmc/")
res_sim <- f_full_run_rep_rjmcmc(list_params_model[[sim_id]])
save(res_sim, file = paste("res_split_id", sim_id, ".RData", sep = ""))
