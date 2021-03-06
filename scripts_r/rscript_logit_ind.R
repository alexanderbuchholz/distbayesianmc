library(distbayesianmc)
#task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#task_id <- as.numeric(task_id_string)
args = commandArgs(trailingOnly=TRUE)
print(args)

i_iter <- as.numeric(args[1])

setwd("~/R_programming/distbayesianmc")
source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
stan_code <- readChar(fileName, file.info(fileName)$size)

set.seed(i_iter)
x_wait <- rexp(1,1)
Sys.sleep(x_wait)

mod <- stan_model(model_code = stan_code, auto_write = T)
#setwd("./sim_results/logistic/")
setwd("/scratch/alexander/distbayesianmc_logit/")
for(dataset in vec_datasets){
  for(ssplits in vec_splits){
    dataset_loaded <- f_dataset_loader(dataset)
    splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=i_iter, typesplit=typesplit)
    splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
    f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = i_iter, iter = i_iter, typesplit = typesplit, nchain = nchain, typeprior=typeprior)
  }
}

