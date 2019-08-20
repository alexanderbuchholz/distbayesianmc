# rscript compile logistic regression model

library(distbayesianmc)
setwd("~/R_programming/distbayesianmc")
source("~/R_programming/distbayesianmc/params_simulation/params_logit_higgs.R")

if (sys.nframe() == 0){
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- rstan::stan_model(model_code = stan_code, auto_write = T, save_dso = FALSE)
  saveRDS(mod, file = "~/R_programming/distbayesianmc/stan_models/fit_logistic.rds")
  rm(mod)
} else {
  mod <- readRDS(file = "~/R_programming/distbayesianmc/stan_models/fit_logistic.rds")
  
  rstan::rstan_options(auto_write = TRUE)
  
  dataset_loaded <- f_dataset_loader("higgs1_small")
  splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=1, iseed=1, typesplit="random")
  splitted_data <- f_prep_prior_logistic(splitted_data, scale = 1)
  res <- f_stan_sampling_single_split(mod, splitted_data[[1]], dataset = "pima", i_seed = 1, iter = 1, typesplit = "random", nchain = 2000, typeprior = "normal")
}
