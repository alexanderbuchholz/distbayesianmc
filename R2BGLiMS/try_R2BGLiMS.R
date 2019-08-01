# R2BGLiMS
library(R2BGLiMS)
library(distbayesianmc)

f_rjmcmc_on_splits <- function(splitted_data_rjmcmc, i_split = 1, n.mil = 2, i_seed=1, thinning.interval=1){
  # bringing the data into the right format
  pvars <- splitted_data_rjmcmc[[i_split]]$Pparams
  X <-  splitted_data_rjmcmc[[i_split]]$X[,2:pvars]
  predictors <- colnames(X)
  data <- cbind(X, splitted_data_rjmcmc[[i_split]]$y)
  colnames(data)[pvars] <- "outcome"
  # preparing the betaprior
  # TODO: make sur that the scale is set right here!
  beta.priors <- cbind(rep(0, pvars-1), rep(splitted_data_rjmcmc[[i_split]]$scale, pvars-1))
  row.names(beta.priors) <- predictors
  
  rjmcmc.results <- R2BGLiMS(
    likelihood="Logistic",
    data=data,
    outcome.var="outcome",
    g.prior = FALSE,
    model.space.priors=list(
      "a"=1, "b"=1, "Variables"=predictors),
    n.mil=n.mil,
    seed=i_seed,
    thinning.interval = thinning.interval,
    beta.priors = beta.priors,
    extra.arguments = list("Adaption_Iterations" = 2e5, 
                           "AlphaPriorSd" = splitted_data_rjmcmc[[i_split]]$scale**0.5
    )
  )
  topmodel_res <- TopModels(rjmcmc.results)  
  mcmc_ouput <- rjmcmc.results@mcmc.output  
  return(list(topmodel_res=topmodel_res, mcmc_ouput=mcmc_ouput, rjmcmc.results=rjmcmc.results))
}
  
prior_prob_k_specific_model <- function(dim_model, dim_all, a=1, b=1){
  # gives the prior probability of a model with a specific dimension, given the full dimension and a and b
  part1 <- beta(dim_model+a, dim_all - dim_model+b)
  part2 <- (dim_all+1)*beta(dim_model+1, dim_all - dim_model+1)*beta(a, b)
  return(part1/(part2 * choose(dim_all,dim_model))) # take into account that we are looking at exactly one single model
}

f_filter_rows_mcmc_rjmcmc_output <- function(mcmc_ouput, topmodel_res, index_model){
  # filter selected colums
  pvars <- dim(topmodel_res)[2]
  selector_model <- c(1, topmodel_res[index_model,(1:pvars-1)])==1
  mcmc_ouput_inter <- mcmc_ouput[,1:pvars]
  filter_matrix <- mcmc_ouput_inter != 0
  selector_rows <- apply(filter_matrix, 1, function(x) all(x==selector_model))
  #selector_rows <- filter_matrix %*% selector_model == sum(selector_mode)
  return(mcmc_ouput_inter[selector_rows, selector_model])
}

f_normal_approx_rjmcmc_output <- function(mcmc_ouput, topmodel_res, index_model){
  # first filter the output
  mcmc_ouput_inter <- f_filter_rows_mcmc_rjmcmc_output(mcmc_ouput, topmodel_res, index_model)
  # calculate means on the filtered ouput
  modelmeans <- colMeans(mcmc_ouput_inter)
  modelcovars <- cov(mcmc_ouput_inter)
  return(list(modelmeans=modelmeans, modelcovars=modelcovars))
}


f_bayesfactor_two_models_rjmcmc <- function(rjmcmcressplit, indexM1, indexM2){
  pvars <- dim(rjmcmcressplit$topmodel_res)[2]
  selectorM1 <- c(1, rjmcmcressplit$topmodel_res[indexM1,(1:pvars-1)])==1
  selectorM2 <- c(1, rjmcmcressplit$topmodel_res[indexM2,(1:pvars-1)])==1
  dimM1 <- sum(selectorM1)-1
  dimM2 <- sum(selectorM2)-1
  
  log_posterior_odds <- log(rjmcmcressplit$topmodel_res[indexM1,pvars]/rjmcmcressplit$topmodel_res[indexM2,pvars])  
  log_bayesfactor <- log_posterior_odds - log(prior_prob_k_specific_model(dimM1, pvars-1, 1, 1)) + log(prior_prob_k_specific_model(dimM2, pvars-1, 1, 1))
  return(log_bayesfactor)
}

f_repeat_rjmcmc_sampling <- function(iterations, splitted_data_rjmcmc){
  # repeat sampling, but just for a single split
  results_topmodels <- list()
  for(i_iter in 1:iterations){
    rjmcmc_split <- f_rjmcmc_on_splits(splitted_data_rjmcmc, i_split = 1, n.mil = 1,i_seed = i_iter, thinning.interval = 100)
    results_topmodels[[i_iter]] <- rjmcmc_split$topmodel_res
  }
  return(results_topmodels)
  
}
scale = 1
ssplits = 1
typesplit = "random"
dataset = "pima"
dataset_loaded <- f_dataset_loader(dataset)
splitted_data_rjmcmc <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=1, typesplit=typesplit)
splitted_data_rjmcmc <- f_prep_prior_logistic(splitted_data_rjmcmc, scale = scale)

results_topmodels <- f_repeat_rjmcmc_sampling(10, splitted_data_rjmcmc)
#rjmcmc_split1 <- f_rjmcmc_on_splits(splitted_data_rjmcmc, i_split = 1, n.mil = 1, thinning.interval = 100)
#rjmcmc_split2 <- f_rjmcmc_on_splits(splitted_data_rjmcmc, i_split = 2, n.mil = 1, thinning.interval = 100)
#TracePlots(rjmcmc_split1$rjmcmc.results, vars.to.include = c("npreg", "glu"))
#AutocorrelationPlot(rjmcmc_split1$rjmcmc.results)

#i_res <-  results_topmodels[[1]]
f_key_for_model_trans_df <- function(i_res_topmodel, counter_sim){
  p_vars <- dim(i_res_topmodel)[2]
  i_res_topmodel <-  data.frame(i_res_topmodel)
  i_res_topmodel$key_model <- apply(i_res_topmodel[,1:p_vars-1],1 , function(x) paste(x, sep="", collapse = ""))
  i_res_topmodel["counter_sim"] <- counter_sim
  return(i_res_topmodel)
}
results_topmodels_backup <- results_topmodels
counter_sim <- 1
for(i_res in results_topmodels){
  i_res <- f_key_for_model_trans_df(i_res, counter_sim)
  results_topmodels[[counter_sim]] <- i_res
  counter_sim <- counter_sim + 1
}

df_sim_res <-  do.call(rbind, results_topmodels)
keys_models <- results_topmodels[[1]]$key_model

indexM1 <- 1
indexM2 <- 4
pvars <- splitted_data_rjmcmc[[1]]$d

#f_bayesfactor_two_models_rjmcmc(rjmcmc_split1, indexM1, indexM2)


selectorM1 <- c(1, results_topmodels_backup[[1]][indexM1,(1:pvars-1)])
selectorM2 <- c(1, results_topmodels_backup[[2]][indexM2,(1:pvars-1)])

source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
stan_code <- readChar(fileName, file.info(fileName)$size)
mod <- stan_model(model_code = stan_code, auto_write = T)

dataset_loaded <- f_dataset_loader(dataset, nobs = nobs)

splitted_dataM1 <- f_pack_split_data(dataset_loaded$X[,selectorM1==1], dataset_loaded$y, ssplits=ssplits, iseed=1, typesplit=typesplit)
splitted_dataM1 <- f_prep_prior_logistic(splitted_dataM1, scale = scale)
res_approx_M1 <-  f_stan_sampling_splitted_data(mod, splitted_dataM1, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
logbf_M1 <- res_approx_M1$normconstcombined
res_approx_M1$vec_logsubpost


splitted_dataM2 <- f_pack_split_data(dataset_loaded$X[,selectorM2==1], dataset_loaded$y, ssplits=ssplits, iseed=1, typesplit=typesplit)
splitted_dataM2 <- f_prep_prior_logistic(splitted_dataM2, scale = scale)
res_approx_M2 <-  f_stan_sampling_splitted_data(mod, splitted_dataM2, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
logbf_M2 <-  res_approx_M2$normconstcombined
res_approx_M2$vec_logsubpost


estimate_single <- logbf_M1 - logbf_M2   + log(prior_prob_k_specific_model(sum(selectorM1)-1, 7, 1, 1)) - log(prior_prob_k_specific_model(sum(selectorM2)-1, 7, 1, 1))

boxplot(log(df_sim_res[df_sim_res$key_model==keys_models[indexM1],]$Post.Prob/df_sim_res[df_sim_res$key_model==keys_models[indexM2],]$Post.Prob))
abline(h=estimate_single, col="red")

res_approx_M1$vec_logsubpost-res_approx_M2$vec_logsubpost

