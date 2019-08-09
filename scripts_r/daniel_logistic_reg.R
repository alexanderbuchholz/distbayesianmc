library(distbayesianmc)

# test one single run

if(T){
  setwd("~/R_programming/distbayesianmc")
  source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
  setwd("./sim_results/logistic/")
  library(doParallel)
  registerDoParallel(cores=1)
  for(dataset in vec_datasets){
    for(ssplits in vec_splits){
      #foreach(iter = 1:iters) %dopar% {
        for(iter in 1:iters){
        
        #dataset_loaded <- f_dataset_loader(dataset)
        #splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=iter, typesplit=typesplit)
        #splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
        f_function_model_split(ssplits, dataset = dataset, iseed = iter, iter = iter, typesplit = typesplit, n_steps = nchain, scale = scale)
        
        
        
      }
    }
  }
}

df <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters)
#f_plot_res_data_frame(df)
f_plot_res_data_frame(df, vec_datasets = vec_datasets)
