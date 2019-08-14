# script evaluate rjmcmc

dataset <-  "higgs1_large"
typesplit <-  "random"
setwd("/scratch/alexander/distbayesianmc_rjmcmc/")
files_sim_results <- list.files()

selector_files <- grepl(dataset, files_sim_results)
list_relevant_files <- files_sim_results[selector_files]

load(list_relevant_files[1])
res_out$mcmc_ouput
