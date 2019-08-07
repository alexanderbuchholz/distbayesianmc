# one split versus two splits

library(distbayesianmc)
library(dplyr)

list_params_model_onesplit <- list(scale =  1,
                          ssplits =  1,
                          typesplit =  "random",
                          dataset =  "pima", 
                          mrep = 18)

list_params_model_multisplits1 <- list(scale =  1,
                                   ssplits =  2,
                                   typesplit =  "strat_y_cluster",
                                   dataset =  "higgs1_small",#"",
                                   mrep = 18)

list_params_model_multisplits2 <- list(scale =  1,
                                       ssplits =  5,
                                       typesplit =  "strat_y_cluster",
                                       dataset =  "higgs1_small",#"",
                                       mrep = 18)

list_params_model_multisplits3 <- list(scale =  1,
                                       ssplits =  10,
                                       typesplit =  "strat_y_cluster",
                                       dataset =  "higgs1_small",#"",
                                       mrep = 18)



res_onesplit <- f_full_run_rep_rjmcmc(list_params_model_onesplit)
save(res_onesplit, file="res_onesplit.RData")

res_multisplit1 <- f_full_run_rep_rjmcmc(list_params_model_multisplits1)
save(res_multisplit1, file="res_multisplit1.RData")

res_multisplit2 <- f_full_run_rep_rjmcmc(list_params_model_multisplits2)
save(res_multisplit2, file="res_multisplit2.RData")

res_multisplit3 <- f_full_run_rep_rjmcmc(list_params_model_multisplits3)
save(res_multisplit3, file="res_multisplit3.RData")




#head(res_onesplit$df_sim_res_all)
#head(res_twosplits$df_sim_res_all)
df1 <- res_onesplit$df_sim_res_all
#df2 <- res_twosplits$df_sim_res_all
keys <- df1[df1$counter_sim==1,]$key_model

index_M1 <- 1
index_M2 <- 2
key1 <- keys[index_M1]
key2 <- keys[index_M2]

bf_rjmcmc <- f_joint_bf_model_splits_rjmcmc(res_onesplit, res_multisplit1, key1, key2)

res_simple_comp <- f_simple_bf_two_models_all_splits(list_params_model_onesplit, res_onesplit$results_splits, index_M1, index_M2)

estimate_single <- res_simple_comp$M1$vec_logsubpost - res_simple_comp$M2$vec_logsubpost + res_simple_comp$M1$model_prior - res_simple_comp$M2$model_prior

pdf(paste("boxplot_m", index_M1, "_m", index_M2, ".pdf" ,sep=""))
boxplot(as.data.frame(bf_rjmcmc))
title(paste("Posterior Odds model ", index_M1, " vs. model ", index_M2, sep=""))
abline(h= estimate_single, col="red")
dev.off()
