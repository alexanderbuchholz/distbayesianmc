# spare linear regression using stan
library(distbayesianmc)




if(T){
  setwd("~/R_programming/distbayesianmc")
  source("~/R_programming/distbayesianmc/params_simulation/params_sparse_median.R")
  stan_code <- readChar(fileName, file.info(fileName)$size)
  
  mod <- rstan::stan_model(model_code = stan_code, auto_write = T)
  setwd("/scratch/alexander/distbayesianmc_sparselinear_median/")
  library(doParallel)
  registerDoParallel(cores=6)
  for(dataset in vec_datasets){
    for(ssplits in vec_splits){
      foreach(i_iter = 1:iters) %dopar% {
        #for(i_iter in 1:20){
        
        dataset_loaded <- f_dataset_loader(dataset, highcorr = highcorr, nobs = nobs)
        splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, Pparams = Pparams, ssplits=ssplits, iseed=i_iter, typesplit=typesplit)
        splitted_data <- f_prep_prior_median(splitted_data, scale = scale)
        f_stan_sampling_splitted_data_median(mod, splitted_data, dataset = dataset, i_seed = i_iter, iter = i_iter, typesplit = typesplit, nchain = nchain, typeprior=typeprior)
        
        
        
      }
    }
  }
}

df <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, 6, type_sim = "stan_median_")
#f_plot_res_data_frame(df)
f_plot_res_data_frame(df, vec_datasets = vec_datasets)
save(df, file="res_sparselinear_ultrasmall_median.Rda")
df %<>% mutate(medianlogsubpost = as.numeric(medianlogsubpost))

p1 <- ggplot(df, aes_string(x="splits", y="medianlogsubpost", fill="dataset")) +
  geom_boxplot() + scale_x_discrete(limits=as.character(sort(as.numeric(levels(df$splits))))) + theme_minimal()
