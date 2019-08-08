# one split versus two splits

library(distbayesianmc)
library(dplyr)

list_params_model_onesplit <- list(scale =  1,
                          ssplits =  1,
                          typesplit =  "random",
                          dataset =  "pima", 
                          mrep = 18,
                          n.mil = 5, 
                          ncores = 4)

list_params_model_multisplits1 <- list(scale =  10,
                                   ssplits =  3,
                                   typesplit =  "random",
                                   dataset =  "pima",#"",
                                   mrep = 12, 
                                   n.mil = 2, 
                                   ncores = 6)



#res_onesplit <- f_full_run_rep_rjmcmc(list_params_model_onesplit)
#save(res_onesplit, file="res_onesplit.RData")

res_multisplit1 <- f_full_run_rep_rjmcmc(list_params_model_multisplits1)
save(res_multisplit1, file="res_multisplit1.RData")




#head(res_onesplit$df_sim_res_all)
#head(res_twosplits$df_sim_res_all)
df1 <- res_multisplit1$df_sim_res_all
#df2 <- res_twosplits$df_sim_res_all
#keys <- df1[df1$counter_sim==1,]$key_model

f_intersection_keys <- function(splits, df){
  keys_per_split1 <- df %>% filter(counter_sim == 1, split==1) %>% dplyr::select(key_model)
  common_keys <- intersect(keys_per_split1, keys_per_split1)    
  for(i_split in 2:splits){
    keys_per_split_i <- df %>% filter(counter_sim == 1, split==i_split) %>% dplyr::select(key_model)
    common_keys <- intersect(common_keys, keys_per_split_i)    
  }
  return(common_keys)
}



common_keys <- f_intersection_keys(list_params_model_multisplits1$ssplits, df1)

i_split <- 2
i_rep <- 3


index_M1 <- 1
index_M2 <- 2
key1 <- common_keys[index_M1,]
key2 <- common_keys[index_M2,]


log(mean( (df1 %>%  filter((key_model == key1  ) & split == i_split) %>% dplyr::select(Post.Prob))[,]  ) / 
  mean((df1 %>%  filter((key_model == key2 ) & split == i_split) %>% dplyr::select(Post.Prob))[,] ) )

#bf_rjmcmc <- f_joint_bf_model_splits_rjmcmc(res_onesplit, res_multisplit1, key1, key2)

res_simple_comp <- f_simple_bf_two_models_all_splits(list_params_model_multisplits1, key1, key2)

estimate_single <- res_simple_comp$M1$vec_logsubpost - res_simple_comp$M2$vec_logsubpost + res_simple_comp$M1$model_prior - res_simple_comp$M2$model_prior

colMeans(res_simple_comp$M1$betasamples[[i_split]])
cov(res_simple_comp$M1$betasamples[[2]])
summary_stats_list1 <- f_summary_stats_per_split(key1, res_multisplit1)
summary_stats_list1[[2]]$mat_means
summary_stats_list1[[2]]$mat_cov

mcmc_output_rjmcmc <- f_filter_rows_mcmc_rjmcmc_output(res_multisplit1$results_splits[[i_split]][[i_rep]]$mcmc_ouput, res_multisplit1$results_splits[[i_split]][[i_rep]]$topmodel_res, key1) 
beta_samples <- res_simple_comp$M1$betasamples[[i_split]]

plot(density(mcmc_output_rjmcmc[,2]))
lines(density(beta_samples[,2]), col="red")


pdf(paste("boxplot_m", index_M1, "_m", index_M2, ".pdf" ,sep=""))
boxplot(as.data.frame(bf_rjmcmc))
title(paste("Posterior Odds model ", index_M1, " vs. model ", index_M2, sep=""))
abline(h= estimate_single, col="red")
dev.off()

