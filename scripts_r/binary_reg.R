library(distbayesianmc)

# test one single run
if(F){
  source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- stan_model(model_code = stan_code, auto_write = T)
  
  dataset_loaded <- f_dataset_loader(dataset, nobs = nobs)
  splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=iter, typesplit=typesplit)
  splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
  res_approx <-  f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = typesplit, nchain = nchain, typeprior=typeprior)
  res_exact <- f_function_model_split(ssplits = ssplits, iter = iter, iseed = iseed, n_steps = nchain, burnin = 200, dataset = dataset, typesplit = typesplit, returnres = T, scale=scale)
  f_plot_grid_params_dens(res_approx)
}

if(F){
  source("~/R_programming/distbayesianmc/params_simulation/params_probit.R")
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- stan_model(model_code = stan_code, auto_write = T)
  
  dataset_loaded <- f_dataset_loader(dataset, nobs = nobs)
  splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=iter, typesplit=typesplit)
  splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
  res_approx <-  f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = typesplit, nchain = nchain, typeprior=typeprior)
  f_plot_grid_params_dens(res_approx)
}


if(T){
  setwd("~/R_programming/distbayesianmc")
  source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
  stan_code <- readChar(fileName, file.info(fileName)$size)
  
  mod <- stan_model(model_code = stan_code, auto_write = T)
  setwd("./sim_results/logistic/")
  library(doParallel)
  registerDoParallel(cores=3)
  for(dataset in vec_datasets){
    for(ssplits in vec_splits){
      #foreach(iter = 1:iters) %dopar% {
        for(iter in 1:iters){
        
        dataset_loaded <- f_dataset_loader(dataset)
        splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=iter, typesplit=typesplit)
        splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
        f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = typesplit, nchain = nchain, typeprior=typeprior)
        
        
        
      }
    }
  }
}

df <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters)
#f_plot_res_data_frame(df)
f_plot_res_data_frame(df, vec_datasets = vec_datasets)
