# spare linear regression using stan
library(distbayesianmc)

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

if(T){
  #setwd("/home/alexander/r_programming/rstan/sim_results_inter")
  setwd("/home/alex/R_programming/rstan/sim_results_inter")
  library(doParallel)
  registerDoParallel(cores=3)
  vec_datasets <- c("sparse_reg_1", "sparse_reg_2")
  #vec_datasets <- c("sparse_reg_1")
  typesplit = "random"
  scale <-  1
  vec_splits <- c(10,20)#, 10,20)#c(2,3,5)#,)#,20,40)
  #vec_splits <- c(20,40)
  iters <- 10
  nobs <- 1000
  for(dataset in vec_datasets){

    for(ssplits in vec_splits){
      foreach(iter = 1:iters) %dopar% {
        #for(iter in 1:20){

        dataset_loaded <- f_dataset_load_logistic_regression(dataset, nobs=nobs, highcorr = T)
        splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, P=2, ssplits=ssplits, iseed=iter, typesplit="random")
        splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)

        res <-  f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = "random", typeprior = "laplace_normal", nchain = 6000)
      }
    }
  }
}


vec_types_splits <- c("random")
df <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters)
f_plot_res_data_frame(df)
f_plot_res_data_frame(df, vec_datasets = vec_datasets)
