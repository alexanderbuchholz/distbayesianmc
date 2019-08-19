# script evaluate rjmcmc

dataset <-  "higgs1_large"
typesplit <-  "random"
setwd("/scratch/alexander/distbayesianmc_rjmcmc/")
#files_sim_results <- list.files()

#selector_files <- grepl(dataset, files_sim_results)
#list_relevant_files <- files_sim_results[selector_files]

#load(list_relevant_files[1])
#res_out$mcmc_ouput

# putting all the model results in one joint frame
vec_splits <- c(1,2,3,5)
dataset <- "sim2"
iters <- 10 # 20
list_all_results <- list()
for (ssplits in vec_splits){
  list_res_splits <- list()
  for(i_split in 1:ssplits){
    list_res_iter <- list()
    for(i_iter in 1:iters){
      sim_res_name <- paste("rjmcmc_output_", dataset, dataset, "_random_splits_", ssplits ,"_i_split_", i_split, "_iter_", i_iter, ".RData", sep = "")
      load(sim_res_name)
      resrjmcmc <- res_out
      list_res_iter[[i_iter]] <- resrjmcmc#$topmodel_res
    }
    df_all_iters <- f_combine_topmodels_in_df(list_res_iter)
    df_all_iters$i_split <- i_split
    df_all_iters$splits <- ssplits
    list_res_splits[[i_split]] <- df_all_iters
  }
  df_sim_res_per_split <-  do.call(rbind, list_res_splits)
  list_all_results[[as.character(ssplits)]] <- df_sim_res_per_split
}

df_all <-  do.call(rbind, list_all_results)
hist(df_all$splits)
