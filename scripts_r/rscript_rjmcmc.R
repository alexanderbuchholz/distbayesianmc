#!/usr/bin/env Rscript
library(distbayesianmc)
#task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#task_id <- as.numeric(task_id_string)
args = commandArgs(trailingOnly=TRUE)
print(args)

sim_id <- as.numeric(args[1])

mrep <- 24
n.mil <- 5
n.cores <- 8
dataset <-  "pima"
typesplit <-  "random"
savestring <-  paste(dataset, typesplit, sep = "_")
idsim <-  0

list_params_model_onesplit <- list(scale =  1,
                                   ssplits =  1,
                                   typesplit =  typesplit,
                                   dataset = dataset, 
                                   mrep = mrep,
                                   n.mil = n.mil, 
                                   n.cores = n.cores, 
                                   savestring = savestring)

list_params_model_multisplits1 <- list(scale =  1,
                                       ssplits =  2,
                                       typesplit =  typesplit,
                                       dataset = dataset,#"",
                                       mrep = mrep, 
                                       n.mil = n.mil,
                                       n.cores = n.cores,
                                       savestring = savestring)

list_params_model_multisplits2 <- list(scale =  1,
                                       ssplits =  3,
                                       typesplit =  typesplit,
                                       dataset = dataset,#"",
                                       mrep = mrep,
                                       n.mil = n.mil,
                                       n.cores = n.cores,
                                       savestring = savestring)

list_params_model_multisplits3 <- list(scale =  1,
                                       ssplits =  5,
                                       typesplit =  typesplit,
                                       dataset = dataset,#"",
                                       mrep = mrep,
                                       n.mil = n.mil,
                                       n.cores = n.cores,
                                       savestring = savestring)

list_params_model_multisplits4 <- list(scale =  1,
                                       ssplits =  10,
                                       typesplit =  typesplit,
                                       dataset = dataset,#"",
                                       mrep = mrep,
                                       n.mil = n.mil,
                                       n.cores = n.cores,
                                       savestring = savestring)


list_params_model <- list(list_params_model_onesplit, 
                          list_params_model_multisplits1, 
                          list_params_model_multisplits2, 
                          list_params_model_multisplits3, 
                          list_params_model_multisplits4)

setwd("/scratch/alexander/distbayesianmc_rjmcmc/")
for(params_model in list_params_model){
  res_sim <- f_single_run_rep_rjmcmc(params_model, sim_id)
  save(res_sim, file = paste("res_split_id", params_model$idsim, "_splits_", params_model$ssplits, ".RData", sep = ""))
}
