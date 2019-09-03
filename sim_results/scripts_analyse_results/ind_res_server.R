# load the distributed data from the server
library(distbayesianmc)

dataset <- "hla1"
ssplits <- 50
i_split <- 1
i_seed <- 1
typesplit <- "random"
file_identifier <- ""
mrep <- 20

setwd("/scratch/alexander/distbayesianmc_sparselinear/res_ind_server/")
#setwd("/scratch/alexander/distbayesianmc_sparselinear/")

vec_datasets <- c("hla1", "hla2")
#vec_datasets <- c("hla_ultra_small1", "hla_ultra_small2")
vec_splits <- c(20, 50, 100)
#vec_splits <- c(1,5,10,20,50)
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


p1 <- ggplot(df , aes_string(x="splits", y="normconstcombined", fill="dataset")) +
  geom_boxplot() +  theme_minimal() +  labs(fill = "Model", y = "log evidence", title="Subset HLA data set") +
  theme(plot.title = element_text(hjust = 0.5, size=18),
        text = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.position="bottom",
  )  +
  scale_x_discrete(labels=c("1 \n(13,235)", "5 \n(2,647)", "10 \n(1,324)", "20 \n(662)","50 \n(265)") )+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5)) +guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  scale_fill_discrete(name = "Model", labels = c("1", "2"))


p2 <- ggplot(df_all , aes_string(x="splits", y="normconstcombined", fill="dataset"))+ #facet_wrap(. ~ dataset, ncol=2)+
  geom_boxplot() +  theme_minimal() +  labs(fill = "Model", y = "", title="Full HLA data set") +
  theme(plot.title = element_text(hjust = 0.5, size=18),
        text = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.position="bottom",
  )  +
  scale_x_discrete(labels=c("20 \n(6,618)", "50 \n(2,647)", "100 \n(1,324)") )+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = c(1.5, 2.5)) +guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  scale_fill_discrete(name = "Model", labels = c("1", "2"))



figure <- ggarrange(p1, p2, widths = c(3.5,2.5),
                    #labels = c("A", "B"),
                    ncol = 2)
ggsave("hladata.pdf", plot = figure, width = 8, height = 4)
