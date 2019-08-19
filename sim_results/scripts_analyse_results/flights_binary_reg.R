# analyse results
library(distbayesianmc)
library(gridExtra)
library(grid)
setwd("~/R_programming/distbayesianmc")
source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
#setwd("./sim_results/logistic/") # use this for the pima dataset
setwd("/scratch/alexander/distbayesianmc_logit/")

df_logistic <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters, type_sim = "stan_")
df_logistic_exact <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits, iters, type_sim = "")
df_logistic_reduced <- df_logistic %>% select(normconstcombined, splits, dataset)
df_logistic_exact_reduced <- df_logistic_exact %>% select(normconstcombined, splits, dataset)


df_logistic_exact_reduced %<>% mutate(model = replace(dataset, dataset=="flights_complex1", "1 exact")) %>% mutate(model = replace(model, model=="flights_complex2", "2 exact"))

df_logistic_reduced %<>% mutate(model = replace(dataset, dataset=="flights_complex1", "1 approx")) %>% mutate(model = replace(model, model=="flights_complex2", "2 approx"))

df_all <- rbind(df_logistic_reduced, df_logistic_exact_reduced)
library(ggplot2)



p1 <- ggplot(df_all , aes_string(x="splits", y="normconstcombined", fill="model")) +
  geom_boxplot() +  theme_minimal() +  labs(fill = "Model", y = "log normalizing constant", title="Estimated normalizing constant \nfor the flights data set") +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        )  + scale_fill_manual(values=c("#0fe600", "#E69F00", "#56B4E9", "#E600B0"))





p2 <- ggplot(df_all  %>% filter(model != "2 exact") %>% filter(splits != "100")  , aes_string(x="splits", y="normconstcombined", fill="model")) +
  geom_boxplot() +  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  ) + scale_fill_manual(values=c("#0fe600", "#E69F00", "#56B4E9"))



g2 <- ggplotGrob(p2)

p3 = p1 + annotation_custom(grob = g2, xmin = 0.5, xmax = 4, ymin=-190000, ymax=-162000)
ggsave("flightsdata.pdf", plot = p3, width = 7, height = 4)
plot(p3)


library(gridExtra)
library(grid)

t1 = arrangeGrob(g1,ncol=1, left = textGrob("A", y = 1, vjust=1, gp=gpar(fontsize=20)))
t2 = arrangeGrob(g2,ncol=1, left = textGrob("B", y = 1, vjust=1, gp=gpar(fontsize=20)))


final = arrangeGrob(t1,t2, layout_matrix = cbind(c(1,2), c(3,3)))
grid.arrange(final)
 


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

