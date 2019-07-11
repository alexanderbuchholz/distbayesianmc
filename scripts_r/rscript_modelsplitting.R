#!/usr/bin/env Rscript
library(distbayesianmc)
#task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#task_id <- as.numeric(task_id_string)
args = commandArgs(trailingOnly=TRUE)
print(args)

ssplits <- as.numeric(args[1])
task_id <- as.numeric(args[2])
model_string <- args[3]
dataset <- args[4]


params_root_string <- "~/R_programming/distbayesianmc/params_simulation/"
source_string <- paste(params_root_string, model_string, sep="")
source(source_string)



# stan model
stan_code <- readChar(fileName, file.info(fileName)$size)
mod <- stan_model(model_code = stan_code, auto_write = T)

setwd("/scratch/alexander/run_sim_modelsplitting/results/")
dataset_loaded <- f_dataset_loader(dataset)
splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed= task_id, typesplit=typesplit)
splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = task_id, iter = task_id, typesplit = typesplit, nchain = nchain, typeprior=typeprior)
