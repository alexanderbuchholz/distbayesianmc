# analyse results
library(distbayesianmc)
library(gridExtra)
library(grid)
setwd("~/R_programming/distbayesianmc")
#source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
#setwd("./sim_results/logistic/") # use this for the pima dataset
setwd("/scratch/alexander/distbayesianmc_logit/")



iters <- 20
vec_datasets <- c("flights_complex1", "flights_complex2")#, )#, "higgs1_full", "higgs2_full")
#vec_datasets <- c("higgs1_large", "higgs2_large")#, "higgs1_full", "higgs2_full")
#vec_datasets <- c("pima")
vec_types_splits_exact <- c("random", "strat_y_cluster")
vec_types_splits_approx <- c("random")

# multiple simulations
#vec_splits <- c(1,2,3,5,10,20)
vec_splits <- c(10,20,50)#,100)

df_logistic <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits_approx, iters, type_sim = "stan_")
df_logistic_exact <- f_combine_const_data_in_frame(vec_splits, vec_datasets, vec_types_splits_exact, iters, type_sim = "")
df_logistic_reduced <- df_logistic %>% select(normconstcombined, splits, dataset, typesplit)
df_logistic_exact_reduced <- df_logistic_exact %>% select(normconstcombined, splits, dataset, typesplit)


df_logistic_exact_reduced %<>% mutate(model = replace(dataset, dataset=="flights_complex1", "1 exact")) %>% mutate(model = replace(model, model=="flights_complex2", "2 exact"))
df_logistic_exact_reduced %<>% mutate(model_strat = replace(typesplit, typesplit=="strat_y_cluster", "stratified")) %>% mutate(model_strat = paste(model_strat, model, sep= " "))

df_logistic_reduced %<>% mutate(model = replace(dataset, dataset=="flights_complex1", "1 approx")) %>% mutate(model = replace(model, model=="flights_complex2", "2 approx"))
df_logistic_reduced %<>% mutate(model_strat = replace(typesplit, typesplit=="strat_y_cluster", "stratified")) %>% mutate(model_strat = paste(model_strat, model, sep= " "))

df_all <- rbind(df_logistic_reduced, df_logistic_exact_reduced)
library(ggplot2)



p1 <- ggplot(df_all %>% filter(splits != "100", typesplit == "random") , aes_string(x="splits", y="normconstcombined", fill="model")) +
  geom_boxplot() +  theme_minimal() +  labs(fill = "", y = "log evidence", title="Estimated evidence flights data, both models") +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        text = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.position="none", # "bottom",
        )  + scale_fill_manual(values=c("#0fe600", "#E69F00", "#56B4E9", "#E600B0")) +
  scale_x_discrete(labels=c("10 \n(32,734)", "20 \n(16,367)", "50 \n(6,547)") )+ 
    theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  + scale_y_continuous(labels = scientific) +
  geom_vline(xintercept = c(1.5, 2.5)) +#guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  annotate("text", x = 0.72, y = -147000, label = "1a", size= 5) + annotate("text", x = 0.89, y = -147000, label = "1e", size= 5) +
  annotate("text", x = 1.07, y = -146500, label = "2a", size= 5) + annotate("text", x = 1.27, y = -147000, label = "2e", size= 5) +
  annotate("text", x = 1.72, y = -147000, label = "1a", size= 5) + annotate("text", x = 1.89, y = -147000, label = "1e", size= 5) +
  annotate("text", x = 2.07, y = -146500, label = "2a", size= 5) + annotate("text", x = 2.27, y = -148000, label = "2e", size= 5) +
  annotate("text", x = 2.72, y = -147000, label = "1a", size= 5) + annotate("text", x = 2.89, y = -147500, label = "1e", size= 5) +
  annotate("text", x = 3.07, y = -146500, label = "2a", size= 5) + annotate("text", x = 3.27, y = -153000, label = "2e", size= 5)

ggsave("flightsdata_comp.pdf", plot = p1, width = 12, height = 6)
p1



p2 <- ggplot(df_all  %>% filter(model != "2 exact", typesplit == "random") %>% filter(splits != "100")  , aes_string(x="splits", y="normconstcombined", fill="model")) +
  geom_boxplot() +  theme_light() + labs(title="Zoom")+
  theme(plot.title = element_text(hjust = 0.5, size=16),
        text = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
  ) + scale_fill_manual(values=c("#0fe600", "#E69F00", "#56B4E9")) +
  scale_x_discrete(labels=c("10 \n(32,734)", "20 \n(16,367)", "50 \n(6,547)") )+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5))

p2

p3 <- ggplot(df_all %>% filter(model != "2 approx", model != "2 exact") , aes_string(x="splits", y="normconstcombined", fill="model_strat")) +
  geom_boxplot() +  theme_minimal() +  labs(fill = "", y = "log evidence", title="Estimated evidence flights data, model 1 only") +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        text = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.position="bottom",
  )  + scale_fill_manual(values=c("#0fe600", "#E69F00", "#56B4E9", "#E600B0")) +
  scale_x_discrete(labels=c("10 \n(32,734)", "20 \n(16,367)", "50 \n(6,547)") )+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5)) +guides(fill=guide_legend(nrow=1,byrow=TRUE))

ggsave("flightsdata_m1.pdf", plot = p3, width = 12, height = 6)
p3
library(ggpubr)
figure <- ggarrange(p1, p3,
                    #labels = c("A", "B"),
                    ncol = 2)
ggsave("flightsdata.pdf", plot = figure, width = 10, height = 7)

g2 <- ggplotGrob(p2)

p3 = p1 + annotation_custom(grob = g2, xmin = 0.5, xmax = 4, ymin=-190000, ymax=-162000)
ggsave("flightsdata.pdf", plot = p3, width = 7, height = 4)
plot(p3)


# calculate normalising constant using IS
dataset_loaded <- f_dataset_loader("flights_complex2")
splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=1, iseed=1, typesplit="random")
splitted_data <- f_prep_prior_logistic(splitted_data, scale = 1)



negf_posterior <- function(x, params.=params) -f_loglik(x, params = params.)-f_logprior(x, params = params.)
resoptim <- optim(splitted_data[[1]]$bprior, negf_posterior, params=splitted_data[[1]], control= c(maxit=10000), hessian=TRUE)


covarproxy <- 2*diag(diag(solve(resoptim$hessian)))
muproxy <- resoptim$par
samples_is <- rmvnorm(n= 200, mean= muproxy, sigma=covarproxy)
logweights <- f_loglik(samples_is, params = splitted_data[[1]])+f_logprior(samples_is, params = splitted_data[[1]]) - dmvnorm(samples_is, mean=muproxy, sigma=covarproxy, log=TRUE)
normconstis <-  log_sum_exp(logweights)-log(1000)
#normconstis <- log(mean(exp(logweights)))
df_all %>% mutate(sqerror = (normconstcombined - normconstis)**2) %>% group_by(dataset, model, splits, model_strat) %>% summarise(MSE = mean(sqerror)**0.5/mean(normconstcombined)*100)



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

