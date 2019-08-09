library(distbayesianmc)
#task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#task_id <- as.numeric(task_id_string)
args = commandArgs(trailingOnly=TRUE)
print(args)

iter <- as.numeric(args[1])

setwd("~/R_programming/distbayesianmc")
source("~/R_programming/distbayesianmc/params_simulation/params_probit.R")
stan_code <- readChar(fileName, file.info(fileName)$size)

mod <- stan_model(model_code = stan_code, auto_write = T)

setwd("/scratch/alexander/distbayesianmc_probit/")
for(dataset in vec_datasets){
  for(ssplits in vec_splits){
    dataset_loaded <- f_dataset_loader(dataset)
    splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=iter, typesplit=typesplit)
    splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
    f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = typesplit, nchain = nchain, typeprior=typeprior)
  }
}

