# spare linear regression using stan
library(distbayesianmc)
if(F){
  ssplits <- 4
  
  iter <- 10
  scale <- 1
  fileName <- "./stan_files/sparse_linear_reg_laplace.stan"
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- stan_model(model_code = stan_code, auto_write = T, verbose=T)
  
  
  
  dataset <- "sparse_reg_2"
  dataset_loaded <- f_dataset_loader(dataset, nobs=4000, highcorr = T)
  splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, P=2, ssplits=ssplits, iseed=iter, typesplit="random")
  splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
  
  res <-  f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = "random", typeprior = "laplace_normal", nchain=4000)
  f_plot_grid_params_dens(res, dataset_loaded$betastar)
  
  # check that the normalising constant is correct
  res$normconstcombined
  index_vars <-  3
  betaval <- res$mat_means[index_vars,]
  cov_diag <- diag(res$mat_cov[,,1])
  n_prop <- 10**5
  var_factor <- 2
  
  log_laplace_prior <- function(betaval, ssplits){
    pparams <- length(betaval)
    res_lik <- sum(dlaplace(betaval[1:(pparams-2)], m = 0, s = scale*ssplits, log = T)) + dnorm(x = betaval[(pparams-1)], mean = 0, sd = (scale*ssplits)**0.5, log = T) +
    dnorm(x = betaval[(pparams)], mean = 0, sd = (scale*ssplits)**0.5, log = T)
    return(res_lik)
  }
  log_lik_sparse <- function(betaval, dataset_loaded){
    pparams <- length(betaval)
    mat_prod <- dataset_loaded$X %*% betaval[1:(pparams-2)]
    y <- dataset_loaded$y
    res_lik <- sum(dnorm(x = (y-mat_prod), mean =  betaval[(pparams-1)], sd = exp(betaval[(pparams)])**0.5, log = T))
    return(res_lik)
  }
  log_posterior <- function(betaval, dataset_loaded, ssplits) log_laplace_prior(betaval, ssplits)+log_lik_sparse(betaval, dataset_loaded)
  #log_posterior(betaval, splitted_data[[index_vars]])
  proposal_samples <- rmvnorm(n = n_prop, mean = betaval, sigma = var_factor*res$mat_cov[,,index_vars])
  log_prop_weights <- dmvnorm(proposal_samples, mean = betaval, sigma = var_factor*res$mat_cov[,,index_vars], log = T)
  log_post_weights <- sapply(1:n_prop, FUN = function(i) log_posterior(proposal_samples[i,], splitted_data[[index_vars]], ssplits))
  log_sum_exp(log_post_weights-log_prop_weights)-log(n_prop)
  res$vec_logsubpost
}




if(T){
  setwd("~/R_programming/distbayesianmc")
  source("~/R_programming/distbayesianmc/params_simulation/params_sparse.R")
  stan_code <- readChar(fileName, file.info(fileName)$size)
  
  mod <- rstan::stan_model(model_code = stan_code, auto_write = T)
  setwd("/scratch/alexander/distbayesianmc_sparselinear/")
  library(doParallel)
  registerDoParallel(cores=6)
  for(dataset in vec_datasets){
    for(ssplits in vec_splits){
      foreach(i_iter = 1:iters) %dopar% {
        #for(i_iter in 1:20){
        
        dataset_loaded <- f_dataset_loader(dataset, highcorr = highcorr, nobs = nobs)
        splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, Pparams = Pparams, ssplits=ssplits, iseed=i_iter, typesplit=typesplit)
        splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
        f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = i_iter, iter = i_iter, typesplit = typesplit, nchain = nchain, typeprior=typeprior)
        
        
        
      }
    }
  }
}

df <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, 10)
#f_plot_res_data_frame(df)
f_plot_res_data_frame(df, vec_datasets = vec_datasets)


if(F){
  df_reduced <- df %>% select(c(normconstcombined, iter, splits, dataset))
  
  
  true_mean <- df_reduced %>% filter(splits == 1)  %>% group_by(splits, dataset) %>% summarise(true_mean = mean(normconstcombined))
  sd_one_split <- df_reduced %>% filter(splits == 1)  %>% group_by(splits, dataset) %>% summarise(std = var(normconstcombined)**0.5)
  
  df_reduced <- df_reduced %>% mutate(true_mean = ifelse(dataset ==  vec_datasets[1], (true_mean %>% filter(dataset ==  vec_datasets[1]))$true_mean, 
                                                         ifelse(dataset ==  vec_datasets[2], (true_mean %>% filter(dataset ==  vec_datasets[2]))$true_mean, NA)))
  
  df_reduced <- df_reduced %>% mutate(sqerror = (true_mean-normconstcombined)**2)
  
  df_true_means <- df_reduced %>% group_by(splits, dataset) %>% summarise(true_mean_all = mean(true_mean))
  
  average_error <- df_reduced %>% group_by(splits, dataset) %>% summarise(MSE = mean(sqerror), VAR = var(normconstcombined))  %>% mutate(srerror = round(MSE**0.5, 3)) 
  
  average_error[['percenterror']] <- round((average_error$srerror/df_true_means$true_mean_all*100), 4)
  
  average_error[["bias_var_ratio"]] <- round(average_error$MSE/average_error$VAR-1, 2)
  
  write.table(t(as.matrix(
    average_error %>% filter(dataset == vec_datasets[1]) %>% select(splits, srerror, percenterror, bias_var_ratio) 
  )), "table1.txt", quote=FALSE, eol="\\\\\n", sep=" & ")
  
  write.table(t(as.matrix(
    average_error %>% filter(dataset == vec_datasets[2]) %>% select(splits, srerror, percenterror, bias_var_ratio) 
  )), "table2.txt", quote=FALSE, eol="\\\\\n", sep=" & ")
  
}
