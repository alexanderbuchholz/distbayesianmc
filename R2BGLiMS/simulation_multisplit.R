# sim multisplit

library(distbayesianmc)
library(dplyr)

list_params_model_onesplit <- list(scale =  1,
                                   ssplits =  1,
                                   typesplit =  "strat_y_cluster",
                                   dataset =  "higgs1_small", 
                                   mrep = 24,
                                   n.mil = 5)

list_params_model_multisplits1 <- list(scale =  1,
                                       ssplits =  2,
                                       typesplit =  "strat_y_cluster",
                                       dataset =  "higgs1_small",#"",
                                       mrep = 24, 
                                       n.mil = 5)

list_params_model_multisplits2 <- list(scale =  1,
                                       ssplits =  5,
                                       typesplit =  "strat_y_cluster",
                                       dataset =  "higgs1_small",#"",
                                       mrep = 24,
                                       n.mil = 5)

list_params_model_multisplits3 <- list(scale =  1,
                                       ssplits =  10,
                                       typesplit =  "strat_y_cluster",
                                       dataset =  "higgs1_small",#"",
                                       mrep = 24,
                                       n.mil = 5)


setwd("/home/alexander/R_programming/distbayesianmc/sim_results/rjmcmc/")
res_onesplit <- f_full_run_rep_rjmcmc(list_params_model_onesplit)
save(res_onesplit, file="res_onesplit.RData")

res_multisplit1 <- f_full_run_rep_rjmcmc(list_params_model_multisplits1)
save(res_multisplit1, file="res_multisplit1.RData")

res_multisplit2 <- f_full_run_rep_rjmcmc(list_params_model_multisplits2)
save(res_multisplit2, file="res_multisplit2.RData")

res_multisplit3 <- f_full_run_rep_rjmcmc(list_params_model_multisplits3)
save(res_multisplit3, file="res_multisplit3.RData")
