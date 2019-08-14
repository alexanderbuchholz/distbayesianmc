# analyse results
library(distbayesianmc)
setwd("~/R_programming/distbayesianmc")
source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
#setwd("./sim_results/logistic/") # use this for the pima dataset
setwd("/scratch/alexander/distbayesianmc_logit/")

df_logistic <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters, type_sim = "stan_")
df_logistic_exact <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters, type_sim = "")

library(ggplot2)
rename_dataset <- function(x){ 
  if(x == "flights_complex1"){
    "1"
  }  
  else {
    "2"
    }
  }
df_logistic %<>% mutate(model = replace(dataset, dataset=="flights_complex1", "1")) %>% mutate(model = replace(model, model=="flights_complex2", "2"))

p1 <- ggplot(df_logistic , aes_string(x="splits", y="normconstcombined", fill="model")) +
  geom_boxplot() +  theme_minimal() +  labs(fill = "Model", y = "log normalizing constant", title="Estimated normalizing constant \nfor the flights data set") +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        )
ggsave("flightsdata.pdf", plot = p1)

# + theme_gray() 


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

