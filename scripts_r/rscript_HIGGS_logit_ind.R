library(distbayesianmc)
library(rstan)
#task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#task_id <- as.numeric(task_id_string)
args = commandArgs(trailingOnly = TRUE)
print(args)
setwd("~/R_programming/distbayesianmc")
source("~/R_programming/distbayesianmc/params_simulation/params_logit_higgs.R")


ssplits <- as.integer(args[2]) # how many splits?

if (args[3] == "split_loop"){
  # split the loop in bits, slurm handles everything
  sim_id <- as.integer(args[1]) # which iteration are we running?
  i_split <- f_split_number(sim_id, ssplits)
  i_iter <- f_i_iter(sim_id, ssplits)
  iters <- 1 # overwrite what comes from the params file
  
} else {
  # looping is done inside R
  i_split <- as.integer(args[1]) # which iteration are we running?
  sim_id <- "loop iters"
}



#set.seed(i_split)
#x_wait <- rexp(1,1)
#Sys.sleep(x_wait)
if (F) {
  # compile the model on every run
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- rstan::stan_model(model_code = stan_code)#, auto_write = T)
}
if (T) {
  # use a precompiled model
  rstan::rstan_options(auto_write = TRUE)
  mod <- readRDS(file = "~/R_programming/distbayesianmc/stan_models/fit_logistic.rds")
}
#setwd("./sim_results/logistic/")
if( Sys.info()[["user"]] == "ab2603" ){
  setwd("/mrc-bsu/scratch/ab2603/distbayesianmc_higgs/")
  server_ind = T
} else {
  setwd("/scratch/alexander/distbayesianmc_higgs/")
  server_ind = F
}


print(tempdir())
#print(mod)

# delete this again! 
#sim_id <- 1
#sim_id <- 498
#sim_id <- 500
#ssplits <-500

#i_split <- f_split_number(sim_id, ssplits)
#i_iter <- f_i_iter(sim_id, ssplits)

# for testing
#sapply(X = c(1:200), FUN = function(x) f_split_number(x, 20))
#sapply(X = c(1:100), FUN = function(x) f_i_iter(x, 20)) == 1

for (dataset in vec_datasets) {
  dataset_loaded <- f_dataset_loader(dataset, server=server_ind)
  for (i_iter in 1:iters) {
    if (args[3] == "split_loop") i_iter <- f_i_iter(sim_id, ssplits)
    splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=i_iter, typesplit=typesplit)
    splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
    print(tempdir())
    print(list.files(tempdir()))
    tryCatch({
      f_stan_sampling_single_split(mod, splitted_data[[i_split]], dataset = dataset, i_seed = i_iter, iter = i_iter, typesplit = typesplit, nchain = nchain, typeprior = typeprior)
      line <-  paste("completed, ", sim_id, ", ", dataset, ",", i_split,",", i_iter, ",", ssplits, sep = "")
      write(line, file = "status_sim.txt", append = TRUE)
      
    }, error = function(err) { 
      line <-  paste("failed, ", sim_id, ", ", dataset, ",", i_split,",", i_iter, ",", ssplits, sep = "")
      write(line, file = "status_sim.txt", append = TRUE)
      
      })
    
  # for (i_iter in 1:iters) {
  #   dataset_loaded <- f_dataset_loader(dataset)
  #   splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=i_iter, typesplit=typesplit)
  #   splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
  #   print(tempdir())
  #   print(list.files(tempdir()))
  #   #tryCatch({
  #   
  #         f_stan_sampling_single_split(mod, splitted_data[[i_split]], dataset = dataset, i_seed = i_iter, iter = i_iter, typesplit = typesplit, nchain = nchain, typeprior = typeprior)
  #       #   },
  #       #   error = function(err) {
  #       #     print("Try exception: recompile model and run again")
  #       #     source("~/R_programming/distbayesianmc/params_simulation/params_logit_higgs.R")
  #       #     print(fileName)
  #       #     stan_code <- readChar(fileName, file.info(fileName)$size)
  #       #     print(stan_code)
  #       #     mod <- rstan::stan_model(model_code = stan_code, save_dso=FALSE)
  #       #     print("compiled, now running")
  #       #     f_stan_sampling_single_split(mod, splitted_data[[i_split]], dataset = dataset, i_seed = i_iter, iter = i_iter, typesplit = typesplit, nchain = nchain, typeprior = typeprior)
  #       #     print("Run finished")
  #       #   }
  #       # )
  #   print(tempdir()) 
  # }
  }
}

