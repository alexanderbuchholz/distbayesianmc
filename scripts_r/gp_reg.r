# gp regression stan

library(rstan)
library(tidyverse)
library(bridgesampling)
source("/home/alexander/r_programming/modelsplitting/R/f_dataset_loader_logistic_reg.R")
source("/home/alexander/r_programming/modelsplitting/R/f_helper.R")
source("/home/alexander/r_programming/modelsplitting/R/functions_logistic_reg.R")
source("/home/alexander/r_programming/modelsplitting/R/functions_product_normal.R")
source("/home/alexander/r_programming/rstan/f_data_splits_stan.R")

fileName <- "./gp_reg.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
mod <- stan_model(model_code = stan_code, auto_write = T)

N = 400
if(F){
  x = 4*pi*runif(N)
  y = cos(x)+rnorm(N)*0.1
  plot(x,y)
  dat <- list(N        = N,
              x        = x,
              y    = y
              )
  
  resStan <- sampling(mod, data = dat, chains = 1, iter = 10000, warmup = 2000, thin = 1, seed = 1)
  bridge_sampler(resStan)
  interres <- rstan::extract(resStan, pars="beta", permuted=F)
  plot(density(interres[,,1]))
  plot(density(interres[,,2]))
  plot(density(interres[,,3]))
}

scale <-  1
dataset <- "gp_test"
ssplits <- 4
iter <- 10
dataset_loaded <- f_dataset_load_logistic_regression(dataset, nobs=N)
splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, d=3, ssplits=ssplits, iseed=iter, typesplit="random")
splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)

res <-  f_stan_sampling_splitted_data_logistic(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = "random", epapprox = F, bridgepack = T)
