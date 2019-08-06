# functions for distributed R2BGLIMS

f_rjmcmc_on_splits <- function(splitted_data_rjmcmc, i_split = 1, n.mil = 2, i_seed=1, thinning.interval=1, savestring=""){
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
      "a"=1, "b"=pvars-1, "Variables"=predictors),
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
  nsplits <-  length(splitted_data_rjmcmc)
  dataset <- splitted_data_rjmcmc[[i_split]]$dataset
  topmodel_res <- f_key_for_model_trans_df(topmodel_res, i_seed)
  
  res_out <- list(topmodel_res=topmodel_res, mcmc_ouput=mcmc_ouput, rjmcmc.results=rjmcmc.results)
  if(savestring != ""){
    fname <- paste("rjmcmc_output_", dataset, "_splits_", nsplits,"_i_split_", i_split,  "_iter_", i_seed, ".RData", sep = "")
    save(res_out, file = fname)
  }
  return(res_out)
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

# function that creates the keys
f_key_for_model_trans_df <- function(i_res_topmodel, counter_sim){
  p_vars <- dim(i_res_topmodel)[2]
  i_res_topmodel <-  data.frame(i_res_topmodel)
  i_res_topmodel$key_model <- apply(i_res_topmodel[,1:p_vars-1],1 , function(x) paste(x, sep="", collapse = ""))
  i_res_topmodel["counter_sim"] <- counter_sim
  return(i_res_topmodel)
}

f_repeat_rjmcmc_sampling <- function(iterations, splitted_data_rjmcmc, i_split = 1, n.mil = 5, savestring = ""){
  # repeat sampling, but just for a single split
  results_topmodels <- list()
  for(i_iter in 1:iterations){
    rjmcmc_split <- f_rjmcmc_on_splits(splitted_data_rjmcmc, i_split = i_split, n.mil = n.mil, i_seed = i_iter, thinning.interval = 100, savestring = savestring)
    #i_res_topmodel <- rjmcmc_split$topmodel_res
    #rjmcmc_split$topmodel_res <- f_key_for_model_trans_df(i_res_topmodel, i_iter)
    results_topmodels[[i_iter]] <- rjmcmc_split
  }
  return(results_topmodels)
  
}

f_parallel_repeat_rjmcmc_sampling <- function(iterations, splitted_data_rjmcmc, i_split = 1, n.mil = 5, ncores = 4, savestring = ""){
  library(doParallel)
  registerDoParallel(cores=ncores)
  
  foreach(iter = 1:iterations) %dopar% {
    #for(iter in 1:20){
    rjmcmc_split <- f_rjmcmc_on_splits(splitted_data_rjmcmc, i_split = i_split, n.mil = n.mil, i_seed = iter, thinning.interval = 100, savestring = savestring)

  }
}





# function that combines the top model results in one single frame
# creates the keys and returns a global frame
f_combine_topmodels_in_df <- function(results_topmodels){
  counter_sim <- 1
  list_res <- list()
  for(i_res in results_topmodels){
    list_res[[counter_sim]] <- i_res$topmodel_res
    counter_sim <- counter_sim + 1
  }
  
  df_sim_res <-  do.call(rbind, list_res)
  return(df_sim_res)
}


f_bf_simple_approx <- function(splitted_data_rjmcmc, results_topmodels, indexM1, indexM2, i_split, typesplit = "random"){
  pvars_all <- splitted_data_rjmcmc[[i_split]]$d  
  pvars <- 1+dim(results_topmodels[[i_split]]$topmodel_res)[2]-3
  selectorM1 <- c(1, as.numeric(results_topmodels[[1]]$topmodel_res[indexM1,(1:(pvars-1))])) # selection of included vars based on the first run of th sampler
  selectorM2 <- c(1, as.numeric(results_topmodels[[1]]$topmodel_res[indexM2,(1:(pvars-1))]))
  
  source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- stan_model(model_code = stan_code, auto_write = T)
  
  
  splitted_dataM1 <- f_pack_split_data(splitted_data_rjmcmc[[i_split]]$X[,selectorM1==1], splitted_data_rjmcmc[[i_split]]$y, ssplits=1, iseed=1, typesplit=typesplit)
  splitted_dataM1 <- f_prep_prior_logistic(splitted_dataM1, scale = splitted_data_rjmcmc[[i_split]]$scale)
  res_approx_M1 <-  f_stan_sampling_splitted_data(mod, splitted_dataM1, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
  logbf_M1 <- res_approx_M1$normconstcombined
  
  splitted_dataM2 <- f_pack_split_data(splitted_data_rjmcmc[[i_split]]$X[,selectorM2==1], splitted_data_rjmcmc[[i_split]]$y, ssplits=1, iseed=1, typesplit=typesplit)
  splitted_dataM2 <- f_prep_prior_logistic(splitted_dataM2, scale = splitted_data_rjmcmc[[i_split]]$scale)
  res_approx_M2 <-  f_stan_sampling_splitted_data(mod, splitted_dataM2, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
  logbf_M2 <-  res_approx_M2$normconstcombined
  estimate_single <- logbf_M1 - logbf_M2   + log(prior_prob_k_specific_model(sum(selectorM1)-1, pvars_all-1, 1, pvars_all-1)) - log(prior_prob_k_specific_model(sum(selectorM2)-1, pvars_all-1, 1, pvars_all-1))  
  return(estimate_single)
}



f_simple_bf_two_models_all_splits <- function(list_params_model, results_topmodels, indexM1, indexM2){
  scale <-  list_params_model$scale
  ssplits <-  list_params_model$ssplits
  typesplit <-  list_params_model$typesplit
  dataset <-  list_params_model$dataset
  dataset_loaded <- f_dataset_loader(dataset)
  
  
  pvars_all <- splitted_data_rjmcmc[[1]]$d  
  pvars <- 1+dim(results_topmodels[[1]]$topmodel_res)[2]-3
  selectorM1 <- c(1, as.numeric(results_topmodels[[1]]$topmodel_res[indexM1,(1:(pvars-1))])) # selection of included vars based on the first run of th sampler
  selectorM2 <- c(1, as.numeric(results_topmodels[[1]]$topmodel_res[indexM2,(1:(pvars-1))]))
  
  source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- stan_model(model_code = stan_code, auto_write = T)
  
  
  splitted_dataM1 <- f_pack_split_data(dataset_loaded$X[,selectorM1==1], dataset_loaded$y, ssplits=ssplits, iseed=1, typesplit=typesplit, dataset = list_params_model$dataset)
  splitted_dataM1 <- f_prep_prior_logistic(splitted_dataM1, scale = list_params_model$scale)
  res_approx_M1 <-  f_stan_sampling_splitted_data(mod, splitted_dataM1, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
  #logbf_M1 <- res_approx_M1$normconstcombined
  res_approx_M1$model_prior <- log(prior_prob_k_specific_model(sum(selectorM1)-1, pvars_all-1, 1, pvars_all-1))
  
  splitted_dataM2 <- f_pack_split_data(dataset_loaded$X[,selectorM2==1], dataset_loaded$y, ssplits=ssplits, iseed=1, typesplit=typesplit, dataset = list_params_model$dataset)
  splitted_dataM2 <- f_prep_prior_logistic(splitted_dataM2, scale = list_params_model$scale)
  res_approx_M2 <-  f_stan_sampling_splitted_data(mod, splitted_dataM2, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
  res_approx_M2$model_prior <- log(prior_prob_k_specific_model(sum(selectorM2)-1, pvars_all-1, 1, pvars_all-1))
  #logbf_M2 <-  res_approx_M2$normconstcombined
  #estimate_single <- logbf_M1 - logbf_M2   + log(prior_prob_k_specific_model(sum(selectorM1)-1, pvars_all-1, 1, pvars_all-1)) - log(prior_prob_k_specific_model(sum(selectorM2)-1, pvars_all-1, 1, pvars_all-1))  
  return(list(M1 = res_approx_M1, M2 = res_approx_M2))
}
