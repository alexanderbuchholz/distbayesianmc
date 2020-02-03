# load the distributed data from the server
library(distbayesianmc)

typesplit <- "random"
file_identifier <- ""
mrep <- 20

setwd("/scratch/alexander/distbayesianmc_higgs/res_ind_server/")

vec_datasets <- c("higgs1_full", "higgs2_full")
vec_splits <- c(100, 500, 1000)
list_splits <- list()
counter <- 1
for(dataset in vec_datasets){
  for(ssplits in vec_splits){
    list_iterations_model <- list()
    for(i_seed in 1:mrep){
      list_results_ind_sim <- list()
      vec_normconsts <- rep(0,ssplits)
      vec_time <- rep(0,ssplits)
      vec_alphasub <- rep(0,ssplits)
      print(paste(dataset, ssplits, i_seed, sep = " "))
      tryCatch({
        
      
        for(i_split in 1:ssplits){
          
          filename_small <- paste("small_sim_stan_ind_", dataset, "_", ssplits, "_splits_i_split_", i_split, "_rep_", i_seed, "_seed_", i_seed, "_", typesplit, file_identifier, ".RData", sep = "") 
          load(file = filename_small)
          if(i_split == 1){
            Pparams <- dim(res_small[["mat_cov"]])[1]
            mat_means <- matrix(0, ncol=Pparams, nrow=ssplits)
            mat_cov <- array(0, dim=c(Pparams, Pparams, ssplits))
            mat_means[i_split,] <- res_small[["mat_means"]]
            mat_cov[,,i_split] <- res_small[["mat_cov"]]  
          }
          else{
            mat_means[i_split,] <- res_small[["mat_means"]]
            mat_cov[,,i_split] <- res_small[["mat_cov"]]  
          }
          vec_normconsts[i_split] <- res_small[["logsubpost"]]
          vec_time[i_split] <- res_small[["runtime"]]
          vec_alphasub[i_split] <- res_small[["alphasub"]]
          
        }
        part_all_normconst <- f_integral_product_gaussian(mat_means, mat_cov)
        list_results_ind_sim[["normconstcombined"]] <- as.numeric(sum(vec_normconsts) + part_all_normconst + mean(vec_alphasub))
        list_results_ind_sim[["varnormconst"]] <- var(vec_normconsts)
        list_results_ind_sim[["meantime"]] <- mean(vec_time)
        list_results_ind_sim[["irep"]] <- i_seed
        list_results_ind_sim[["splits"]] <- ssplits
        list_results_ind_sim[["dataset"]] <- dataset
      }, error = function(cond){
        list_results_ind_sim[["normconstcombined"]] <- NA
        list_results_ind_sim[["varnormconst"]] <- NA
        list_results_ind_sim[["meantime"]] <- NA
        list_results_ind_sim[["irep"]] <- i_seed
        list_results_ind_sim[["splits"]] <- ssplits
        list_results_ind_sim[["dataset"]] <- dataset
        
      })
      df_ind_res <- as.data.frame(list_results_ind_sim)
      list_iterations_model[[i_seed]] <- df_ind_res
    }
    
    list_splits[[counter]] <- do.call(rbind, list_iterations_model)
    counter <- counter + 1 
  }
}


df_all <- do.call(rbind, list_splits)

df_all$splits <- as.factor(df_all$splits)
#p1 <- ggplot(df_all, aes_string(x="splits", y="normconstcombined", fill="dataset")) +
#  geom_boxplot() + scale_x_discrete(limits=as.character(sort(as.numeric(levels(df_all$splits))))) + theme_minimal()# + theme_gray()


df_all %<>% mutate(model = recode(dataset, higgs1_full = "1",
                                 higgs2_full = "2"))


p1 <- ggplot(df_all , aes_string(x="splits", y="normconstcombined", fill="model")) +
  geom_boxplot() +  theme_minimal() +  labs(fill = "Model", y = "log evidence", title="Estimated evidence \nfor the Higgs data set") +
  theme(plot.title = element_text(hjust = 0.5, size=24),
        text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        legend.position="bottom",
  )  +
  scale_x_discrete(labels=c("100 \n(110k)", "500 \n(22k)", "1000 \n(11k)") )+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5)) +guides(fill=guide_legend(nrow=1,byrow=TRUE))+ scale_y_continuous(labels = scientific)
#+ scale_colour_Publication()+ theme_Publication()



p2 <- ggplot(df_all , aes_string(x="splits", y="normconstcombined") )+ facet_wrap(. ~ model, ncol=2, scales="free_y")+
  geom_boxplot() +  theme_minimal() +  labs(fill = "Model", y = "log evidence", title="Estimated evidence") +
  theme(plot.title = element_text(hjust = 0.5, size=24),
        text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        legend.position="bottom",
  )  +
  scale_x_discrete(labels=c("100 \n(110k)", "500 \n(22k)", "1000 \n(11k)") )+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5)) +guides(fill=guide_legend(nrow=2,byrow=TRUE)) #+ scale_y_continuous(labels = scientific)#+scale_colour_Publication()+ theme_Publication()


library(ggpubr)
figure_combined <- ggarrange(p1, p2, nrow=1)

ggsave("/scratch/alexander/higgsdata_combined.pdf", plot = figure_combined, width = 15, height = 5)


ggsave("higgsdata.pdf", plot = p2, width = 8, height = 4)

p3 <- ggplot(df_all , aes_string(x="splits", y="meantime") )+ facet_wrap(. ~ model, ncol=1, scales="free_y")+
  geom_boxplot() +  theme_minimal() +  labs(fill = "Model", y = "Time (secs)", title="Average computation time per shard") +
  theme(plot.title = element_text(hjust = 0.5, size=24),
        text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        legend.position="bottom",
  )  +
  scale_x_discrete(labels=c("100 \n(110,000)", "500 \n(22,000)", "1000 \n(11,000)") )+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5)) +guides(fill=guide_legend(nrow=2,byrow=TRUE))

