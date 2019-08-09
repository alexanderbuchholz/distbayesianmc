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
  
  topmodel_res <- TopModels(rjmcmc.results, n.top.models = 100, remove.empty.cols = F)  
  mcmc_ouput <- rjmcmc.results@mcmc.output  
  nsplits <-  length(splitted_data_rjmcmc)
  dataset <- splitted_data_rjmcmc[[i_split]]$dataset
  topmodel_res <- f_key_for_model_trans_df(topmodel_res, i_seed)
  
  res_out <- list(topmodel_res=topmodel_res, mcmc_ouput=mcmc_ouput, rjmcmc.results=rjmcmc.results)
  if(savestring != ""){
    fname <- paste("rjmcmc_output_", dataset, savestring, "_splits_", nsplits,"_i_split_", i_split,  "_iter_", i_seed, ".RData", sep = "")
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
  # check if the beta binomial prior is correct!
  log_posterior_odds <- log(rjmcmcressplit$topmodel_res[indexM1,pvars]/rjmcmcressplit$topmodel_res[indexM2,pvars])  
  log_bayesfactor <- log_posterior_odds - log(prior_prob_k_specific_model(dimM1, pvars-1, 1, pvars-1)) + log(prior_prob_k_specific_model(dimM2, pvars-1, 1, pvars-1))
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

f_repeat_rjmcmc_sampling <- function(iterations, splitted_data_rjmcmc, i_split = 1, n.mil = 5, savestring = "", i_seed = 1){
  # repeat sampling, but just for a single split
  results_topmodels <- list()
  for(i_iter in 1:iterations){
    rjmcmc_split <- f_rjmcmc_on_splits(splitted_data_rjmcmc, i_split = i_split, n.mil = n.mil, i_seed = i_seed, thinning.interval = 100, savestring = savestring)
    #i_res_topmodel <- rjmcmc_split$topmodel_res
    #rjmcmc_split$topmodel_res <- f_key_for_model_trans_df(i_res_topmodel, i_iter)
    results_topmodels[[i_iter]] <- rjmcmc_split
  }
  return(results_topmodels)
  
}

f_parallel_repeat_rjmcmc_sampling <- function(iterations, splitted_data_rjmcmc, i_split = 1, n.mil = 5, ncores = 4, savestring = ""){
  library(doParallel)
  library(parallel)
  #cl <- parallel::makeCluster(ncores)
  #doParallel::registerDoParallel(cl)
  registerDoParallel(cores=ncores)
  
  foreach(iter = 1:iterations) %dopar% {
    #for(iter in 1:20){
    #rjmcmc_split <- f_rjmcmc_on_splits(splitted_data_rjmcmc, i_split = i_split, n.mil = n.mil, i_seed = iter, thinning.interval = 100, savestring = savestring)
    distbayesianmc::f_rjmcmc_on_splits(splitted_data_rjmcmc, i_split = i_split, n.mil = n.mil, i_seed = iter, thinning.interval = 100, savestring = savestring)

  }
  #parallel::stopCluster(cl)
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

# OLD FUNCTION DO NOT USE ANYMORE

# f_bf_simple_approx <- function(splitted_data_rjmcmc, results_topmodels, indexM1, indexM2, i_split, typesplit = "random"){
#   pvars_all <- splitted_data_rjmcmc[[i_split]]$d  
#   pvars <- 1+dim(results_topmodels[[i_split]]$topmodel_res)[2]-3
#   selectorM1 <- c(1, as.numeric(results_topmodels[[1]]$topmodel_res[indexM1,(1:(pvars-1))])) # selection of included vars based on the first run of th sampler
#   selectorM2 <- c(1, as.numeric(results_topmodels[[1]]$topmodel_res[indexM2,(1:(pvars-1))]))
#   
#   source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
#   stan_code <- readChar(fileName, file.info(fileName)$size)
#   mod <- stan_model(model_code = stan_code, auto_write = T)
#   
#   
#   splitted_dataM1 <- f_pack_split_data(splitted_data_rjmcmc[[i_split]]$X[,selectorM1==1], splitted_data_rjmcmc[[i_split]]$y, ssplits=1, iseed=1, typesplit=typesplit)
#   splitted_dataM1 <- f_prep_prior_logistic(splitted_dataM1, scale = splitted_data_rjmcmc[[i_split]]$scale)
#   res_approx_M1 <-  f_stan_sampling_splitted_data(mod, splitted_dataM1, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
#   logbf_M1 <- res_approx_M1$normconstcombined
#   
#   splitted_dataM2 <- f_pack_split_data(splitted_data_rjmcmc[[i_split]]$X[,selectorM2==1], splitted_data_rjmcmc[[i_split]]$y, ssplits=1, iseed=1, typesplit=typesplit)
#   splitted_dataM2 <- f_prep_prior_logistic(splitted_dataM2, scale = splitted_data_rjmcmc[[i_split]]$scale)
#   res_approx_M2 <-  f_stan_sampling_splitted_data(mod, splitted_dataM2, dataset = dataset, i_seed = 1, iter = 1, typesplit = typesplit, nchain = 10000, typeprior="normal")
#   logbf_M2 <-  res_approx_M2$normconstcombined
#   estimate_single <- logbf_M1 - logbf_M2   + log(prior_prob_k_specific_model(sum(selectorM1)-1, pvars_all-1, 1, pvars_all-1)) - log(prior_prob_k_specific_model(sum(selectorM2)-1, pvars_all-1, 1, pvars_all-1))  
#   return(estimate_single)
# }



f_simple_bf_two_models_all_splits <- function(list_params_model, key1, key2){
  #browser()
  ssplits <-  list_params_model$ssplits
  typesplit <-  list_params_model$typesplit
  dataset <-  list_params_model$dataset
  dataset_loaded <- f_dataset_loader(dataset)
  
  
  pvars_all <- dim(dataset_loaded$X)[2] #splitted_data_rjmcmc[[1]]$d  
  #pvars <- 1+dim(results_topmodels[[1]][[1]]$topmodel_res)[2]-3
  
  selectorM1 <- c(1, as.numeric(unlist(strsplit(key1, "")))) # selection of included vars based on the first run of th sampler
  selectorM2 <- c(1, as.numeric(unlist(strsplit(key2, ""))))
  
  source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- rstan::stan_model(model_code = stan_code, auto_write = T)
  #browser()
  
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

f_full_run_rep_rjmcmc <- function(list_params_model){
  dataset_loaded <- f_dataset_loader(list_params_model$dataset)
  splitted_data_rjmcmc <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=list_params_model$ssplits, iseed=1, typesplit=list_params_model$typesplit, dataset = list_params_model$dataset)
  splitted_data_rjmcmc <- f_prep_prior_logistic(splitted_data_rjmcmc, scale = list_params_model$scale)
  #i_split <-  1
  mrep <- list_params_model$mrep
  # repeat the sampling on one single split
  results_splits <- list()
  results_splits_topmodel_frame <- list()
  for(i_split in 1:list_params_model$ssplits){
    if(iflist_params_model$n.cores==1){
      results_topmodels_parallel <- f_repeat_rjmcmc_sampling(mrep, splitted_data_rjmcmc, i_split = i_split, n.mil = list_params_model$n.mil, savestring = list_params_model$savestring)
    }
    else{
      results_topmodels_parallel <- f_parallel_repeat_rjmcmc_sampling(mrep, splitted_data_rjmcmc, i_split = i_split, n.mil = list_params_model$n.mil, ncores = list_params_model$n.cores, savestring = list_params_model$savestring)
    }
    results_splits[[i_split]] <- results_topmodels_parallel
    df_sim_res <- f_combine_topmodels_in_df(results_topmodels_parallel)
    df_sim_res[["split"]] <- i_split
    results_splits_topmodel_frame[[i_split]] <- df_sim_res
  }
  df_sim_res_all <-  do.call(rbind, results_splits_topmodel_frame)
  return(list(results_splits=results_splits, df_sim_res_all = df_sim_res_all))
}

f_single_run_rep_rjmcmc <- function(list_params_model, i_seed){
  dataset_loaded <- f_dataset_loader(list_params_model$dataset)
  splitted_data_rjmcmc <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=list_params_model$ssplits, iseed=1, typesplit=list_params_model$typesplit, dataset = list_params_model$dataset)
  splitted_data_rjmcmc <- f_prep_prior_logistic(splitted_data_rjmcmc, scale = list_params_model$scale)
  #i_split <-  1
  
  # repeat the sampling on one single split
  results_splits <- list()
  results_splits_topmodel_frame <- list()
  for(i_split in 1:list_params_model$ssplits){
    results_topmodels_parallel <- f_repeat_rjmcmc_sampling(1, splitted_data_rjmcmc, i_split = i_split, n.mil = list_params_model$n.mil, savestring = list_params_model$savestring, i_seed = i_seed)
    #results_splits[[i_split]] <- results_topmodels_parallel
    #df_sim_res <- f_combine_topmodels_in_df(results_topmodels_parallel)
    #df_sim_res[["split"]] <- i_split
    #results_splits_topmodel_frame[[i_split]] <- df_sim_res
  }
  #df_sim_res_all <-  do.call(rbind, results_splits_topmodel_frame)
  #return(list(results_splits=results_splits, df_sim_res_all = df_sim_res_all))
}


f_filter_rows_mcmc_rjmcmc_output <- function(mcmc_output, topmodel_res, key_model){
  # filter selected colums
  pvars <- dim(topmodel_res)[2]-2 # remove key model and counter sim
  selector_model <- c(1, topmodel_res[topmodel_res$key_model == key_model,(1:pvars-1)])==1
  mcmc_ouput_inter <- mcmc_output[,1:pvars]
  filter_matrix <- mcmc_ouput_inter != 0
  selector_rows <- apply(filter_matrix, 1, function(x) all(x==selector_model))
  #selector_rows <- filter_matrix %*% selector_model == sum(selector_mode)
  return(mcmc_ouput_inter[selector_rows, selector_model])
}

f_summary_stats_per_split <- function(key, res_split_list){
  results_splits <- res_split_list$results_splits
  splits <- length(results_splits)
  mrep <- length(results_splits[[1]])
  summary_stats_list <- list()
  pvars <- sum(unlist(strsplit(key, "")) == "1") + 1 # number of variables
  #browser()
  
  summary_stats_list <-  list()
  for(i_rep in 1:mrep){
    # loop over repetitions
    mat_means <- matrix(0, ncol=pvars, nrow=splits)
    mat_cov <- array(0, dim=c(pvars, pvars, splits))
    summary_stats_list[[i_rep]] <- list()
    # loop over the splits
    for(i_split in 1:splits){
      mcmc_output_inter <- f_filter_rows_mcmc_rjmcmc_output(results_splits[[i_split]][[i_rep]]$mcmc_ouput, results_splits[[i_split]][[i_rep]]$topmodel_res, key) 
      mat_means[i_split,] <- colMeans(mcmc_output_inter)
      mat_cov[,,i_split] <- cov(mcmc_output_inter)
    }
    summary_stats_list[[i_rep]]$mat_means <- mat_means
    summary_stats_list[[i_rep]]$mat_cov <- mat_cov
  }
  return(summary_stats_list) # first layer: splits, second layer iterations
}



f_joint_bf_model_splits_rjmcmc <- function(res_onesplit, res_several_splits, key1, key2, list_params_model_multisplits){
  mrep <- length(res_onesplit$results_splits[[1]])
  bfonesplit <- rep(0, mrep)
  
  # extract means and variance of splits running over repetions
  summary_stats_list1 <- f_summary_stats_per_split(key1, res_several_splits)
  summary_stats_list2 <- f_summary_stats_per_split(key2, res_several_splits)
  
  df1 <- res_onesplit$df_sim_res_all
  df2 <- res_several_splits$df_sim_res_all
  
  bfseveralsplits <- rep(0, mrep)
  browser()
  for(i_rep in 1:mrep){
    # single split
    #browser()
    postprob1 <- df1 %>% filter(key_model == key1 | key_model == key2) %>% filter(counter_sim == i_rep) %>% dplyr::select("Post.Prob") 
    bfonesplit[i_rep] <- log(postprob1[1,]/postprob1[2,])
    
    # several splits
    M1dim <- sum(unlist(strsplit(key1, "")) == "1") + 1 # number of variables
    M2dim <- sum(unlist(strsplit(key2, "")) == "1") + 1 # number of variables
    fulldim_no_intercept <- nchar(key1)
    
    df2 %>% filter(key_model == key1 | key_model == key2) %>% filter(counter_sim == i_rep) %>% dplyr::select(Post.Prob, split, key_model)
    df2_key1 <- df2 %>% filter(key_model == key1) %>% filter(counter_sim == 1) %>% dplyr::select(Post.Prob, split)
    df2_key2 <- df2 %>% filter(key_model == key2) %>% filter(counter_sim == 1) %>% dplyr::select(Post.Prob, split)
    bfsplits <- df2_key1$Post.Prob/df2_key2$Post.Prob
    
    log_alpha_key1 <- list_params_model_multisplits$ssplits*f_alpha_sub(ssplits = list_params_model_multisplits$ssplits, Vprior = diag(rep(list_params_model_multisplits$scale, M1dim)), typeprior = "normal")
    log_alpha_key2 <- list_params_model_multisplits$ssplits*f_alpha_sub(ssplits = list_params_model_multisplits$ssplits, Vprior = diag(rep(list_params_model_multisplits$scale, M2dim)), typeprior = "normal")
    
    log_priorM1 <- log(prior_prob_k_specific_model(M1dim-1, fulldim_no_intercept, 1, fulldim_no_intercept))
    log_priorM2 <- log(prior_prob_k_specific_model(M2dim-1, fulldim_no_intercept, 1, fulldim_no_intercept))
    log_prior_odds <- log_priorM1 - log_priorM2
    
    prod1 <- f_integral_product_gaussian(summary_stats_list1[[i_rep]]$mat_means, summary_stats_list1[[i_rep]]$mat_cov)
    prod2 <- f_integral_product_gaussian(summary_stats_list2[[i_rep]]$mat_means, summary_stats_list2[[i_rep]]$mat_cov)
    
    
    log_bf_splits <- sum(log(bfsplits)) - (list_params_model_multisplits$ssplits-1)*log_prior_odds + log_alpha_key1 - log_alpha_key2 + prod1 - prod2
    #log(bffull)
    
    bfseveralsplits[i_rep] <- log_bf_splits
  }
  return(list(bfonesplit=bfonesplit, bfseveralsplits=bfseveralsplits))
  
  
}

