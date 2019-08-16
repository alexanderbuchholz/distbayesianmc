library(distbayesianmc)
#task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#task_id <- as.numeric(task_id_string)
args = commandArgs(trailingOnly = TRUE)
print(args)

i_split <- as.integer(args[1]) # which iteration are we running?
ssplits <- as.integer(args[2]) # how many splits?

setwd("~/R_programming/distbayesianmc")
source("~/R_programming/distbayesianmc/params_simulation/params_logit_higgs.R")
stan_code <- readChar(fileName, file.info(fileName)$size)

#set.seed(i_split)
#x_wait <- rexp(1,1)
#Sys.sleep(x_wait)

mod <- stan_model(model_code = stan_code)#, auto_write = T)
#setwd("./sim_results/logistic/")
setwd("/scratch/alexander/distbayesianmc_higgs/")
print(tempdir())
#print(mod)
for (dataset in vec_datasets) {
  for (i_iter in 1:iters) {
    dataset_loaded <- f_dataset_loader(dataset)
    splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=i_iter, typesplit=typesplit)
    splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
    tryCatch({
          f_stan_sampling_single_split(mod, splitted_data[[i_split]], dataset = dataset, i_seed = i_iter, iter = i_iter, typesplit = typesplit, nchain = nchain, typeprior = typeprior)
          },
          error = function(err) {
            print("Try exception: recompile model and run again")
            source("~/R_programming/distbayesianmc/params_simulation/params_logit_higgs.R")
            print(fileName)
            stan_code <- readChar(fileName, file.info(fileName)$size)
            print(stan_code)
            mod <- stan_model(model_code = stan_code, save_dso=FALSE)
            print("compiled, now running")
            f_stan_sampling_single_split(mod, splitted_data[[i_split]], dataset = dataset, i_seed = i_iter, iter = i_iter, typesplit = typesplit, nchain = nchain, typeprior = typeprior)
            print("Run finished")
          }
        )
    print(tempdir()) 
  }
}

