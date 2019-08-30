# script evaluate rjmcmc

library(distbayesianmc)

typesplit <-  "random"
setwd("/scratch/alexander/distbayesianmc_rjmcmc/")
source("~/R_programming/distbayesianmc/scripts_r/rscript_rjmcmc.R")

#setwd("~/R_programming/exchange_files_server/rjmcmc_sim2/")
#files_sim_results <- list.files()

#selector_files <- grepl(dataset, files_sim_results)
#list_relevant_files <- files_sim_results[selector_files]

#load(list_relevant_files[1])
#res_out$mcmc_ouput

# putting all the model results in one joint frame
vec_splits <- c(1,2,3,5)#,10)
dataset <- "sim2"
iters <-  20

# use this only for the higgs data set
list_params_model <- list_params_model[1:4]

list_all_results_df <- list()
list_all_results_output <- list()
list_all_results_all <- list()
for (ssplits in vec_splits){
  list_res_splits_df <- list()
  list_res_splits_output <- list()
  for(i_split in 1:ssplits){
    list_res_iter <- list()
    for(i_iter in 1:iters){
      sim_res_name <- paste("rjmcmc_output_", dataset, dataset, "_random_splits_", ssplits ,"_i_split_", i_split, "_iter_", i_iter, ".RData", sep = "")
      load(sim_res_name)
      resrjmcmc <- res_out
      list_res_iter[[i_iter]] <- resrjmcmc #$topmodel_res
    }
    df_all_iters <- f_combine_topmodels_in_df(list_res_iter)
    df_all_iters$i_split <- i_split
    df_all_iters$splits <- ssplits
    list_res_splits_df[[i_split]] <- df_all_iters
    list_res_splits_output[[i_split]] <- list_res_iter
  }
  df_sim_res_per_split <-  do.call(rbind, list_res_splits_df)
  list_all_results_all[[as.character(ssplits)]] <- list(results_splits=list_res_splits_output, df_sim_res_all=df_sim_res_per_split)
  list_all_results_df[[as.character(ssplits)]] <- df_sim_res_per_split
  list_all_results_output[[as.character(ssplits)]] <- list_res_splits_output
}

# list_all_results_output 
# nested list: first level: type of splits, second level: individual splits, third level: interations
df_all <-  do.call(rbind, list_all_results_df)

common_keys_single <- f_intersection_keys(list_params_model_onesplit$ssplits, list_all_results_df[['1']])#, df_single_split)

for (params_model in list_params_model){
  common_keys_multisplit <- f_intersection_keys(params_model$ssplits, list_all_results_df[[as.character(params_model$ssplits)]])#, df_single_split)
  keys_intersect <- intersect(common_keys_single, common_keys_multisplit)
}

combos <- list(c(1,2))#, c(2,3), c(3,4), c(1,3), c(2,4), c(1,4))
list_res <- list()
counter_combs <- 0
for(model_comb in combos){
  counter_combs <- counter_combs + 1
  key1 <- keys_intersect[model_comb[1],]
  key2 <- keys_intersect[model_comb[2],]
  
  counter <-  1
  frames_list <- list()
  for (params_model in list_params_model){
    if(counter != 1){
      bf_rjmcmc <- f_joint_bf_model_splits_rjmcmc(list_all_results_all$`1`, list_all_results_all[[as.character(params_model$ssplits)]], key1, key2, params_model)
      df_bf_rjmcmc <- as_tibble(bf_rjmcmc)
      df_bf_rjmcmc %<>% select(bfseveralsplits)
      colnames(df_bf_rjmcmc)[1] <- "BF"
      df_bf_rjmcmc %<>% mutate(splits = params_model$ssplits) %>% mutate(model = paste(model_comb, collapse = "/"))
      
      #boxplot(df_bf_rjmcmc)
      frames_list[[counter-1]] <- df_bf_rjmcmc
    }
    counter <-  counter + 1
  }
  # handle the single split
  df_bf_rjmcmc <- as_tibble(bf_rjmcmc)
  df_bf_rjmcmc %<>% select(bfonesplit)
  colnames(df_bf_rjmcmc)[1] <- "BF"
  df_bf_rjmcmc %<>% mutate(splits = 1) %>% mutate(model = paste(model_comb, collapse = "/"))
  frames_list[[counter-1]] <- df_bf_rjmcmc
  frames_combined <- do.call(rbind, frames_list)
  list_res[[counter_combs]] <- frames_combined
}

frames_different_models <- do.call(rbind, list_res)
library(ggplot2)
frames_different_models %<>% mutate(ssplits = as.character(splits))

frames_different_models$ssplits <- factor(frames_different_models$ssplits, levels = vec_splits, ordered = TRUE)
bp <- ggplot(frames_different_models, aes(y=BF, x=ssplits, group=splits)) + 
  geom_boxplot() + facet_wrap(. ~ model, ncol=3) +  theme_minimal() +  
  labs(fill = "Model", y = "log Bayes factor", x = "# splits", title="log Bayes factors accross \n different splits and models") +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
  ) 

#bp

ggsave(filename = "rjmcmc.pdf", plot = bp, width = 7, height = 4)
