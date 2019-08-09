# analyse results

setwd("~/R_programming/distbayesianmc")
source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
setwd("./sim_results/logistic/")

df_logistic <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters, type_sim = "stan_")
df_logistic_exact <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters, type_sim = "")


setwd("~/R_programming/distbayesianmc")
source("~/R_programming/distbayesianmc/params_simulation/params_probit.R")
setwd("./sim_results/probit/")

df_probit <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters)

library(dplyr)
df_logistic$model <- "logit"
df_probit$model <- "probit"
df_logistic_exact$model <- "logit_ex"
df_logistic_reduced <- df_logistic %>% select(normconstcombined, model, splits, iter)
df_probit_reduced <-df_probit %>% select(normconstcombined, model, splits, iter)
df_logistic_exact <- df_logistic_exact %>% select(normconstcombined, model, splits, iter)

df_all <- rbind(df_logistic_reduced, df_probit_reduced, df_logistic_exact)

df_all %<>% mutate(splitsnew = as.character(splits))

library(ggplot2)
p1 <- ggplot(df_all %>% filter(splitsnew != "20"), aes_string(x="splits", y="normconstcombined", fill="model")) +
  geom_boxplot() + scale_x_discrete(limits=as.character(sort(as.numeric(levels(df$splits))))) + theme_minimal()# + theme_gray()

