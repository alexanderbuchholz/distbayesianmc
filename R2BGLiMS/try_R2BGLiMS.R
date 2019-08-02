# R2BGLiMS
#library(R2BGLiMS)
library(distbayesianmc)

scale = 1
ssplits = 20
typesplit = "random"
dataset = "sim2"
dataset_loaded <- f_dataset_loader(dataset)
splitted_data_rjmcmc <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=1, typesplit=typesplit)
splitted_data_rjmcmc <- f_prep_prior_logistic(splitted_data_rjmcmc, scale = scale)
i_split = 2
# repeat the sampling on one single split
results_topmodels <- f_repeat_rjmcmc_sampling(5, splitted_data_rjmcmc, i_split = i_split, n.mil = 2)
#TracePlots(rjmcmc_split1$rjmcmc.results, vars.to.include = c("npreg", "glu"))
#AutocorrelationPlot(rjmcmc_split1$rjmcmc.results)

df_sim_res <- f_combine_topmodels_in_df(results_topmodels)

# the keys for the different models
keys_models <- df_sim_res[df_sim_res$counter_sim==1,]$key_model

indexM1 <- 2
indexM2 <- 3


#f_bayesfactor_two_models_rjmcmc(rjmcmc_split1, indexM1, indexM2)

f_bf_simple_approx <- function(splitted_data_rjmcmc, results_topmodels, indexM1, indexM2, i_split){
  pvars <- splitted_data_rjmcmc[[i_split]]$d  
  selectorM1 <- c(1, as.numeric(results_topmodels[[1]][indexM1,(1:pvars-1)])) # selection of included vars based on the first run of th sampler
  selectorM2 <- c(1, as.numeric(results_topmodels[[1]][indexM2,(1:pvars-1)]))
  
  source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- stan_model(model_code = stan_code, auto_write = T)
  
  typesplit = "random"
  splitted_dataM1 <- f_pack_split_data(splitted_data_rjmcmc[[i_split]]$X[,selectorM1==1], splitted_data_rjmcmc[[i_split]]$y, ssplits=1, iseed=1, typesplit=typesplit)
  splitted_dataM1 <- f_prep_prior_logistic(splitted_dataM1, scale = splitted_data_rjmcmc[[i_split]]$scale)
  res_approx_M1 <-  f_stan_sampling_splitted_data(mod, splitted_dataM1, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
  logbf_M1 <- res_approx_M1$normconstcombined
  
  splitted_dataM2 <- f_pack_split_data(splitted_data_rjmcmc[[i_split]]$X[,selectorM2==1], splitted_data_rjmcmc[[i_split]]$y, ssplits=1, iseed=1, typesplit=typesplit)
  splitted_dataM2 <- f_prep_prior_logistic(splitted_dataM2, scale = splitted_data_rjmcmc[[i_split]]$scale)
  res_approx_M2 <-  f_stan_sampling_splitted_data(mod, splitted_dataM2, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
  logbf_M2 <-  res_approx_M2$normconstcombined
  estimate_single <- logbf_M1 - logbf_M2   + log(prior_prob_k_specific_model(sum(selectorM1)-1, pvars-1, 1, 1)) - log(prior_prob_k_specific_model(sum(selectorM2)-1, pvars-1, 1, 1))  
  return(estimate_single)
}



estimate_single <- f_bf_simple_approx(splitted_data_rjmcmc, results_topmodels, indexM1, indexM2, i_split)


boxplot(log(df_sim_res[df_sim_res$key_model==keys_models[indexM1],]$Post.Prob/df_sim_res[df_sim_res$key_model==keys_models[indexM2],]$Post.Prob))
abline(h=estimate_single, col="red")

