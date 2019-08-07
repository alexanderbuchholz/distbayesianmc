# R2BGLiMS
#library(R2BGLiMS)
library(distbayesianmc)

list_params_model <- list(scale =  1,
                          ssplits =  1,
                          typesplit =  "random",
                          dataset =  "higgs1_small")

dataset_loaded <- f_dataset_loader(list_params_model$dataset)
splitted_data_rjmcmc <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=list_params_model$ssplits, iseed=1, typesplit=list_params_model$typesplit, dataset = list_params_model$dataset_loaded$dataset)
splitted_data_rjmcmc <- f_prep_prior_logistic(splitted_data_rjmcmc, scale = list_params_model$scale)
#i_split <-  1
mrep <- 10
# repeat the sampling on one single split
results_splits <- list()
results_splits_topmodel_frame <- list()
for(i_split in 1:list_params_model$ssplits){
  results_topmodels_parallel <- f_parallel_repeat_rjmcmc_sampling(mrep, splitted_data_rjmcmc, i_split = i_split, n.mil = 2, ncores = 6, savestring = "test")
  results_splits[[i_split]] <- results_topmodels_parallel
  df_sim_res <- f_combine_topmodels_in_df(results_topmodels_parallel)
  df_sim_res[["split"]] <- i_split
  results_splits_topmodel_frame[[i_split]] <- df_sim_res
}

# TODO: 
# 1. merge the frames for all splits based on RJMCMC
df_sim_res_all <-  do.call(rbind, results_splits_topmodel_frame)
# 2. filter rows that have been most visited
keys_models <- df_sim_res_all[df_sim_res_all$counter_sim==1,]$key_model
indexM1 <- 1
indexM2 <- 3

key1 <- keys_models[indexM1]
key2 <- keys_models[indexM2]
df_sim_res_all_M12 <- df_sim_res_all[(df_sim_res_all$key == key1 | df_sim_res_all$key == key2),]

# 3. make plot using ggplot2
# calculate the ratio first

# make a grouped boxplot


#results_topmodels <- f_repeat_rjmcmc_sampling(2, splitted_data_rjmcmc, i_split = i_split, n.mil = 1)
#TracePlots(rjmcmc_split1$rjmcmc.results, vars.to.include = c("npreg", "glu"))
#AutocorrelationPlot(rjmcmc_split1$rjmcmc.results)

#df_sim_res <- f_combine_topmodels_in_df(results_topmodels_parallel)
#colnames_sim_res <- colnames(df_sim_res)
# the keys for the different models
#keys_models <- df_sim_res[df_sim_res$counter_sim==1,]$key_model


res_simple_comp <- f_simple_bf_two_models_all_splits(list_params_model, results_topmodels_parallel, indexM1, indexM2)
estimate_single <- res_simple_comp$M1$vec_logsubpost - res_simple_comp$M2$vec_logsubpost + res_simple_comp$M1$model_prior - res_simple_comp$M2$model_prior
#f_bayesfactor_two_models_rjmcmc(rjmcmc_split1, indexM1, indexM2)


#estimate_single <- f_bf_simple_approx(splitted_data_rjmcmc, results_topmodels_parallel, indexM1, indexM2, i_split, typesplit = typesplit)

spec_split <- 1
boxplot(log(df_sim_res_all_M12[df_sim_res_all_M12$key_model==keys_models[indexM1] & df_sim_res_all_M12$split==spec_split,]$Post.Prob/df_sim_res_all_M12[df_sim_res_all_M12$key_model==keys_models[indexM2] & df_sim_res_all_M12$split==spec_split,]$Post.Prob))
abline(h=estimate_single, col="red")


#f_normal_approx_rjmcmc_output(results_splits[[1]][[1]]$mcmc_ouput, topmodel_res = results_splits_topmodel_frame, index_model = 1)

