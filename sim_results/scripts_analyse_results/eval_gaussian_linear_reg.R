# evaluate gaussian linear reg

library(distbayesianmc)
library(gridExtra)
library(grid)
setwd("~/R_programming/distbayesianmc")
setwd("/scratch/alexander/distbayesianmc_sparselinear/")



iters <- 20
vec_datasets <- c("gaussian1", "gaussian2", "gaussian3","gaussian4", "gaussian5", "gaussian6") # 
vec_types_splits <- c(typesplit)

# multiple simulations
vec_splits <- c(1,2,5,10,20,50)#1,30,50)#,10,20,50)
typesplit <-  "random"

# multiple simulations
load(file="res_gaussianlinear.Rda")
library(ggplot2)



p1 <- ggplot(df, aes_string(x="splits", y="normconstcombined", fill="dataset")) +
  geom_boxplot() +  theme_minimal() +  labs(fill = "dataset", y = "log evidence", title="Estimated evidence synthetic gaussian model") +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        text = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.position="bottom", # "bottom",
  )  +
  scale_x_discrete(labels=c("1 \n(10,000)", "2 \n(5,000)", "5 \n(2,000)", "10 \n(1,000)", "20 \n(500)", "50 \n(200)") )+ 
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + geom_vline(xintercept = c(1.5, 2.5, 3.5,4.5,5.5)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  scale_fill_discrete(name = "Model", labels = c("1", "2", "3", "4", "5", "6")) 
  
p1


# do a bayesfactor plot
#head(df)
res_list <- list()
i <- 1
for(split in vec_splits){
  df_inter <- df %>% filter(splits == split)
  df_bf_inter <- df_inter %>%
    group_by(dataset) %>%
    mutate(BF = normconstcombined - filter(df_inter, dataset == 'gaussian6')$normconstcombined)
  res_list[[i]] <- df_bf_inter
  i <- i+1
}
df_bf <- bind_rows(res_list)
df_bf$BF <-  df_bf$BF*-1

p2 <- ggplot(df_bf %>% filter(dataset != 'gaussian6'), aes_string(x="splits", y="BF", fill="dataset")) +
  geom_boxplot() +  theme_minimal() +  labs(fill = "dataset", y = "log BF", title="Log BF synthetic Gaussian model: 6 vs. rest") +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        text = element_text(size=18),
        axis.title.x = element_text(size=18) ,# element_blank(),# ,
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=18),
        legend.position = "top"
  )  +
  scale_x_discrete(labels=c("1 \n(10,000)", "2 \n(5,000)", "5 \n(2,000)", "10 \n(1,000)", "20 \n(500)", "50 \n(200)") )+ 
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  + geom_vline(xintercept = c(1.5, 2.5, 3.5,4.5,5.5)) +
  scale_fill_discrete(name = "BF", labels = c("6/1", "6/2", "6/3", "6/4", "6/5")) +
guides(fill=guide_legend(nrow=1,byrow=TRUE))

p2


# contributions approximate term
df %<>% mutate(logsubposteriorapprox = as.numeric(logsubposteriorapprox))
res_list <- list()
i <- 1
for(split in vec_splits){
  df_inter <- df %>% filter(splits == split)
  df_cont_inter <- df_inter %>%
    group_by(dataset) %>%
    mutate(contribution = -logsubposteriorapprox/normconstcombined *100)
  res_list[[i]] <- df_cont_inter
  i <- i+1
}
df_cont <- bind_rows(res_list)


p3 <- ggplot(df_cont %>% filter(splits != "1") , aes_string(x="splits", y="contribution", fill="dataset")) +
  geom_boxplot() +  theme_minimal() +  labs(fill = "dataset", y = "Contribution in %", title = TeX("Contribution of $\\int \\prod_{s=1}^S \\tilde{p}(\\theta | y_s) d \\theta$ to total estimated evidence")) +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        text = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.position="bottom", # "bottom",
  )  +
  scale_x_discrete(labels=c("2 \n(5,000)", "5 \n(2,000)", "10 \n(1,000)", "20 \n(500)", "50 \n(200)") )+ #"1 \n(10,000)", 
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + geom_vline(xintercept = c(1.5, 2.5, 3.5,4.5,5.5)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  scale_fill_discrete(name = "Model", labels = c("1", "2", "3", "4", "5", "6")) 


p3


# plot introduced error by log approx subposterior


if(T){
  df_reduced <- df %>% select(c(normconstcombined, iter, splits, dataset))
  df_reduced %<>% group_by(splits, dataset) %>% mutate(true_mean = mean(normconstcombined))
  df_reduced %<>% mutate(sqerror = (normconstcombined-true_mean)**2)
  df_summary <- df_reduced  %>% group_by(splits, dataset) %>%  summarise(pmse = -(mean(sqerror)**0.5)/mean(true_mean) * 100)
  df_summary %<>%  ungroup() %>% mutate(splits = as.numeric(as.character(splits)))
  df_summary %<>% mutate(obs = 10000/splits)
}

#df_summary %<>%  ungroup() %>% mutate(splits = as.numeric(splits))
p4 <- ggplot(df_summary ) +
  geom_line(aes_string(x="splits", y="pmse", color = "dataset"), size=1.2) +  theme_minimal() +  labs(fill = "dataset", y = "Percentage MSE", title = "Introduced bias from approximation") +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        text = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.position="bottom", # "bottom",
  )  +
  #scale_x_discrete(labels=c("2 \n(5,000)", "5 \n(2,000)", "10 \n(1,000)", "20 \n(500)", "50 \n(200)") )+ #"1 \n(10,000)", 
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + #+ geom_vline(xintercept = c(1.5, 2.5, 3.5,4.5,5.5)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))  +
  #scale_fill_discrete(name = "Model", labels = c("1", "2", "3", "4", "5", "6")) 
  scale_color_discrete(name = "Model", labels = c("1", "2", "3", "4", "5", "6"))

p4

library(ggpubr)
figure1 <- ggarrange( p2, p1,
                    #labels = c("A", "B"),
                    ncol = 1)
figure2 <- ggarrange( p3, p4,
                      #labels = c("A", "B"),
                      ncol = 2)

ggsave("BF_gaussian.pdf", plot = figure1, width = 15, height = 10)
ggsave("approx_error_gaussian.pdf", plot = figure2, width = 15, height = 7)



