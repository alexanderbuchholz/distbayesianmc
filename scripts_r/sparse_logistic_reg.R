library(distbayesianmc)


fileName <- "./stan_files/sparse_logistic_reg_laplace.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
mod <- stan_model(model_code = stan_code, auto_write = T)

if(F){
  scale <-  1
  dataset <- "pima"
  ssplits <- 5
  iter <- 42
  dataset_loaded <- f_dataset_loader(dataset, nobs = 10**3)
  splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=iter, typesplit="random")
  splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
  
  res <-  f_stan_sampling_splitted_data_logistic(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = "random", epapprox = F, bridgepack = T, typeprior = "laplace", nchain = 2000)
  f_plot_grid_params_dens(res)
}

if(T){
  setwd("/home/alexander/r_programming/rstan/sim_results")
  library(doParallel)
  registerDoParallel(cores=6)
  dataset = "sim1"
  typesplit = "random"
  scale <-  1
  vec_splits <- c(1,2,3,5,10)#,10,20)
  
  for(ssplits in vec_splits){
    foreach(iter = 1:20) %dopar% {
      #for(iter in 1:20){
      
      dataset_loaded <- f_dataset_load_logistic_regression(dataset)
      splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=iter, typesplit=typesplit)
      splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
      f_stan_sampling_splitted_data_logistic(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = typesplit)
      
      
      
    }
  }
}
