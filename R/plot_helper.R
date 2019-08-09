# functions for plotting of the results

f_plot_grid_params_dens <- function(res, betastar=NA){
  ssplits <- length(res$betasamples)
  frame_list <- list()
  for(i_s in 1:ssplits){
    df_samples_inter <- data.frame(res$betasamples[[i_s]])
    df_samples_inter$split <- as.character(i_s)
    frame_list[[i_s]] <- df_samples_inter
    #plot(density(res$betasamples[[i_s]][,1]))
  }
  
  beta_frame <-  do.call(rbind, frame_list)
  col_names <- names(beta_frame)
  num_coefs <- length(col_names)-1
  nCol <- floor(sqrt(num_coefs))
  plist <- list()
  for(i_coef in 1:num_coefs){
    p_inter <- ggplot(beta_frame, aes_string(x=col_names[i_coef], fill=col_names[num_coefs+1])) + 
      geom_density(alpha=0.3) + theme(legend.position = "none") 
    if(!is.na(betastar)) p_inter <-  p_inter + geom_vline(xintercept = betastar[i_coef])
    plist[[i_coef]] <- p_inter
  }
  p <- do.call("grid.arrange", c(plist, ncol=nCol))
  p
}


f_combine_const_data_in_frame <- function(vec_splits, vec_datasets, vec_types_splits, iters, type_sim = "stan_"){
  vec_iter <- c(1:iters)
  data_list <- list()
  icounter <- 1
  for(dataset in vec_datasets){
    for(ssplits in vec_splits){
      for(iter in vec_iter){
        for(typesplit in vec_types_splits){
          filename_small <- paste("small_sim_", type_sim, dataset, "_", ssplits, "_splits_rep_", iter, "_seed_" , iter, "_", typesplit, ".RData",sep="")
          #browser()
          load(file=filename_small)
          res_small[["iter"]] <- iter
          res_small[["splits"]] <- ssplits
          res_small[["typesplit"]] <- typesplit
          res_small[["var_logsubposteriors"]] <- var(res_small$vec_logsubpost)
          res_small[["vec_logsubpost"]] <- NULL
          res_small[["vec_loggaussianconst"]] <- NULL
          res_small[["vec_logsubpost_ep"]] <- NULL
          res_small[["normconstcombined"]] <- as.numeric(res_small[["normconstcombined"]])
          
          res_small[["mat_cov"]] <- NULL
          res_small[["mat_means"]] <- NULL
          res_small[["vec_condintegral"]] <- NULL
          res_small$vec_logsubposteriors <- NULL
          res_small$vec_logsubpost_is <- NULL
          res_small$betasamples <- NULL
          res_small[["dataset"]] <- dataset
          #res_to_keep <- list(iter=iter, ssplits=ssplits, normconstcombinedapprox=as.vector(res_small[["normconstcombinedapprox"]]))
          #data_list[[icounter]] <- res_to_keep
          data_list[[icounter]] <- res_small
          icounter <- icounter+1
        }
      }
    }
  }
  #browser()
  df <- as_tibble(matrix(unlist(data_list), nrow=length(data_list), byrow=T))
  names(df) <- names(res_small)
  f_convert_factor <- function(x) as.numeric(levels(x))[x]
  f_convert_to_numeric <- function(x) as.numeric(gsub("," ,".", x))
  df %<>% mutate(splits = as_factor(splits))
  
  df %<>% mutate(var_logsubposteriors =  f_convert_to_numeric(var_logsubposteriors))
  df %<>% mutate(normconstcombined =  f_convert_to_numeric(normconstcombined))
  df
}

# change the normconstcombined back !
f_plot_res_data_frame <- function(df, vec_datasets = NA, vec_types_splits = NA){
  f_name <- paste("normconst", paste(vec_datasets,  collapse = ''), ".pdf", sep="")
  if(all(!is.na(vec_datasets))){
  p1 <- ggplot(df, aes_string(x="splits", y="normconstcombined", fill="dataset")) +
    geom_boxplot() + scale_x_discrete(limits=as.character(sort(as.numeric(levels(df$splits))))) + theme_minimal()# + theme_gray()
  }
  else if(!is.na(vec_types_splits)){
    p1 <- ggplot(df, aes_string(x="splits", y="normconstcombined", fill="typesplit")) +
      geom_boxplot() + scale_x_discrete(limits=as.character(sort(as.numeric(levels(df$splits))))) + theme_minimal()# + theme_gray()
  }
  else {
    p1 <- ggplot(df, aes_string(x="splits", y="normconstcombined")) +
      geom_boxplot() + scale_x_discrete(limits=as.character(sort(as.numeric(levels(df$splits))))) + theme_minimal()# + theme_gray()
  }
  weird <- scales::trans_new("signed_log",
                             transform=function(x) sign(x)*log(abs(x)),
                             inverse=function(x) sign(x)*exp(abs(x)))
  p1 <-  p1 + theme(text = element_text(size=20)) +
    scale_y_continuous(trans=weird)
  
  ggsave(f_name, plot=p1, device = pdf())
}
