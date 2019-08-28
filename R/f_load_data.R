# load data sets

f_sim_sparse_data <- function(nobs, highcorr=F){
  list_data <- list()
  set.seed(123)
  N <-  nobs
  P <-  50
  
  alpha <-  3
  sigma <-  1
  sig_prob <-  0.2
  
  
  x <- matrix(0, ncol = P, nrow = N)
  y <- rep(0, N)
  betavec <- rep(0, P)
  
  # simulate betavec
  for (m in 1:(P-1)) {
    if (rbinom(1, prob = sig_prob, size=1)){
      if (rbinom(1, prob = 0.5, size=1)){
        betavec[m] <-  rnorm(n=1, mean=10, sd=1)
      }
      else{
        betavec[m] <-  rnorm(n=1, mean=-10, sd=1)
      }
    }
    else{
      betavec[m] <-  rnorm(n=1, mean=0, sd=0.25)
    }
  }
  betavec[P] <-  rnorm(n=1, mean=10, sd=1)
  # simulate data
  if(highcorr) {
    corr_mat <- matrix(0.99, nrow = P, ncol = P)
    diag(corr_mat) <- rep(1,P)
    
    var_x <- rep(1,P)
    sigma_x <- diag(var_x) %*% corr_mat %*% diag(var_x)
  }
  else sigma_x <- diag(nrow = P)
  for (n in 1:N) {
    x[n,] <- mvtnorm::rmvnorm(1, sigma=sigma_x)
    #for (m in 1:P){
    #  x[n, m] = rnorm(n=1)
    #}
    y[n] = rnorm(n=1, mean = t(x[n,]) %*% betavec + alpha, sd = sigma)
  }
  
  #list_data[["X"]] <- x[,1:(P-1)] # remove the last observation here
  #list_data[["y"]] <- y
  #list_data[["betastar"]] <- c(beta[1:(P-1)], alpha, log(sigma))
  list_data[["X"]] <- x # remove the last observation here
  list_data[["y"]] <- y
  list_data[["betastar"]] <- c(betavec, alpha, log(sigma))
  list_data[["dataset"]] <- dataset
  return(list_data)
  
}
f_dataset_loader <- function(dataset="pima", nobs=5*10**3, highcorr = T, server=F){
  list_data <- list()
  if(dataset == "pima"){
    data_pima = rbind(Pima.tr, Pima.te)
    X = scale(data_pima[,1:7])
    #X = data_pima[,1:7]
    ones = rep(1, nrow(data_pima))
    X = cbind(ones, X)
    y = (data_pima[,8] == "Yes")*1
    list_data[["X"]] <- X
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "gp_test"){
    set.seed(123)
    N = nobs
    x = 4*pi*runif(N)
    y = cos(x)+rnorm(N)*0.1
    
    
    list_data[["X"]] <- matrix(x, nrow = N, ncol = 1) # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "sparse_reg_1"){
    list_data <- f_sim_sparse_data(nobs, highcorr)
    P <- dim(list_data$X)[2]
    list_data$X <- list_data$X[,1:(P-1)]
    selector <- rep(T, (P + 2))
    selector[P] <- F
    list_data[["betastar"]] <- list_data[["betastar"]][selector]
    list_data[["dataset"]] <- dataset
    #browser()
  }
  else if (dataset == "sparse_reg_2"){
    list_data <- f_sim_sparse_data(nobs, highcorr)
    #list_data[["dataset"]] <- dataset
    
  }
  else if (dataset == "hla"){
    if( F ) {
      # preprocessing the data, only to run once
      df_genotypes <- read.csv("/project/wtccc/ILIKE/hla_common_additive0.8.tsv/hla_common_additive0.8.tsv", header = T, sep = "")
      df_genotypes <- df_genotypes %>% select(-c("NA.", "NA..1"))
      df_outcome <- read.csv("/project/wtccc/ILIKE/hla_for_daniel/ukbb_ukbil_merged.sample", header = T, sep = "", stringsAsFactors = F)
      write.csv(df_genotypes, file = "/scratch/alexander/hladata/hla_genotypes_full.csv")
      
      df_outcome = df_outcome[-1,]
      df_outcome <- df_outcome %>% mutate(mcv_gwas_normalised = as.numeric(mcv_gwas_normalised))
      #df_genotypes <- df_genotypes %>% select(-c("NA.", "NA..1"))
      df_genotypes[['mcv_gwas_normalised']] <- df_outcome$mcv_gwas_normalised
      write.csv(df_genotypes, file = "/scratch/alexander/hladata/hla_genotypes_full_outcome.csv")
      
      df_genotypes <- df_genotypes %>% sample_frac(0.1)
      
      write.csv(df_genotypes, file = "/scratch/alexander/hladata/hla_genotypes_frac_outcome.csv")
      
      df_genotypes <- df_genotypes %>% drop_na()
      ncol <- dim(df_genotypes)[2]
      selector_rows <- (runif(ncol-1) < 0.1)
      df_genotypes_small <- df_genotypes[,which(c(selector_rows, T),T)]
      # corr_vec <- rep(0,ncol-1)
      # for(i in 1:(ncol-1)){
      #   corr_vec[i] <- cor(df_genotypes[,i], df_genotypes[,ncol])
      # }
      write.csv(df_genotypes_small, file = "/scratch/alexander/hladata/hla_genotypes_subset_frac_outcome.csv")
      genenames <- colnames(df_genotypes_small)[1:(dim(df_genotypes_small)[2]-1)]
      reg1 <- lm(as.formula(paste("mcv_gwas_normalised~", paste(genenames, collapse="+"))), df_genotypes_small)
      library(glmnet)
      
  
    }
    
    if(F){
    df_full <- read.csv("/scratch/alexander/hladata/hla_genotypes_full_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    df_full <- df_full %>% dplyr::select(-c("X")) %>% drop_na()
    
    #y <- df_full$mcv_gwas_normalised
    #X <- as.matrix(df_full %>% dplyr::select(-c("mcv_gwas_normalised")))

    ncol <- dim(df_full)[2]-1
    nrow <- dim(df_full)[1]
    
    corr_vec <- rep(0,ncol)
    for(i in 1:(ncol)){
      corr_vec[i] <- cor(df_full[,i], df_full$mcv_gwas_normalised)
    }
    sort_res <- sort(log(abs(corr_vec)), decreasing = T, index.return=T)
    X_reduced1 <- df_full[,sort_res$ix[1:50]]
    X_reduced1 <- cbind(rep(1,nrow), X_reduced1)
    colnames(X_reduced1)[1] <- "intercept"
    X_reduced1[["mcv_gwas_normalised"]] <- df_full$mcv_gwas_normalised
    write.csv(X_reduced1, file = "/scratch/alexander/hladata/hla_1.csv")
    
    X_reduced2 <- df_full[,sort_res$ix[1:100]]
    X_reduced2 <- cbind(rep(1,nrow), X_reduced2)
    colnames(X_reduced2)[1] <- "intercept"
    X_reduced2[["mcv_gwas_normalised"]] <- df_full$mcv_gwas_normalised
    write.csv(X_reduced2, file = "/scratch/alexander/hladata/hla_2.csv")
    
    list_data[["X"]] <- X_reduced # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
    }
    
  }
  else if (dataset == "hla1"){
    if(server){
      df_small <- read.csv("/mrc-bsu/scratch/ab2603/hladata/hla_1.csv", header = T, sep = ",", stringsAsFactors = F)
    }
    else{
      df_small <- read.csv("/scratch/alexander/hladata/hla_1.csv", header = T, sep = ",", stringsAsFactors = F)
    }
    
    #df_small <- read.csv("~/R_programming/exchange_files_server/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    df_small <- df_small %>% dplyr::select(-c("X"))
    y <- df_small$mcv_gwas_normalised
    X <- as.matrix(df_small %>% dplyr::select(-c("mcv_gwas_normalised")))
    
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "hla2"){
    if(server){
      df_small <- read.csv("/mrc-bsu/scratch/ab2603/hladata/hla_2.csv", header = T, sep = ",", stringsAsFactors = F)
    }
    else{
      df_small <- read.csv("/scratch/alexander/hladata/hla_2.csv", header = T, sep = ",", stringsAsFactors = F)
    }
    
    #df_small <- read.csv("~/R_programming/exchange_files_server/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    df_small <- df_small %>% dplyr::select(-c("X"))
    y <- df_small$mcv_gwas_normalised
    X <- as.matrix(df_small %>% dplyr::select(-c("mcv_gwas_normalised")))
    
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "hla_small"){
    df_small <- read.csv("/scratch/alexander/hladata/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    #df_small <- read.csv("~/R_programming/exchange_files_server/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    df_small <- df_small %>% dplyr::select(-c("X"))
    y <- df_small$mcv_gwas_normalised
    X <- as.matrix(df_small %>% dplyr::select(-c("mcv_gwas_normalised")))
    ncol <- dim(X)[2]
    nrow <- dim(X)[1]
    corr_vec <- rep(0,ncol-1)
    for(i in 1:(ncol)){
       corr_vec[i] <- cor(X[,i], y)
    }
    sort_res <- sort(log(abs(corr_vec)), decreasing = T, index.return=T)
    X_reduced <- X[,sort_res$ix[1:50]]
    X_reduced <- cbind(rep(1,nrow), X_reduced)
    colnames(X_reduced)[1] <- "intercept"
    list_data[["X"]] <- X_reduced # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "hla_ultra_small1"){
    if(server){
      df_small <- read.csv("/mrc-bsu/scratch/ab2603/hladata/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    }
    else{
      df_small <- read.csv("/scratch/alexander/hladata/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    }
    
    #df_small <- read.csv("~/R_programming/exchange_files_server/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    df_small <- df_small %>% dplyr::select(-c("X"))
    nsamples <- 10000
    y <- df_small$mcv_gwas_normalised[1:nsamples]
    X <- as.matrix(df_small %>% dplyr::select(-c("mcv_gwas_normalised")))
    X <- X[1:nsamples,]
    ncol <- dim(X)[2]
    nrow <- dim(X)[1]
    corr_vec <- rep(0,ncol)
    for(i in 1:(ncol)){
      corr_vec[i] <- cor(X[,i], y)
    }
    sort_res <- sort(log(abs(corr_vec)), decreasing = T, index.return=T)
    X_reduced <- X[,sort_res$ix[1:20]]
    X_reduced <- cbind(rep(1,nrow), X_reduced)
    colnames(X_reduced)[1] <- "intercept"
    list_data[["X"]] <- X_reduced # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "hla_ultra_small2"){
    if(server){
      df_small <- read.csv("/mrc-bsu/scratch/ab2603/hladata/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    }
    else{
      df_small <- read.csv("/scratch/alexander/hladata/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    }
    
    #df_small <- read.csv("~/R_programming/exchange_files_server/hla_genotypes_subset_frac_outcome.csv", header = T, sep = ",", stringsAsFactors = F)
    df_small <- df_small %>% dplyr::select(-c("X"))
    nsamples <- 10000
    y <- df_small$mcv_gwas_normalised[1:nsamples]
    X <- as.matrix(df_small %>% dplyr::select(-c("mcv_gwas_normalised")))
    X <- X[1:nsamples,]
    ncol <- dim(X)[2]
    nrow <- dim(X)[1]
    corr_vec <- rep(0,ncol)
    for(i in 1:(ncol)){
      corr_vec[i] <- cor(X[,i], y)
    }
    sort_res <- sort(log(abs(corr_vec)), decreasing = T, index.return=T)
    X_reduced <- X[,sort_res$ix[1:40]]
    X_reduced <- cbind(rep(1,nrow), X_reduced)
    colnames(X_reduced)[1] <- "intercept"
    list_data[["X"]] <- X_reduced # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs1"){
    # subset of the higgs data set
    df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1.csv", header = F)
    samplesize = dim(df)[1]
    X = cbind(rep(1, samplesize), df[,23:29])
    colnames(X)[1] <- "V1"
    y = df[,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs1_full"){
    # subset of the higgs data set
    if(server){
      df <- read.csv("/mrc-bsu/scratch/ab2603/higgsdata/HIGGS.csv", header = F)
    }
    else{
      df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv", header = F)
    }
    
    samplesize = dim(df)[1]
    X = cbind(rep(1, samplesize), df[,23:29])
    colnames(X)[1] <- "V1"
    y = df[,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs1_small"){
    # subset of the higgs data set
    if(server){
      df <- read.csv("/mrc-bsu/scratch/ab2603/higgsdata/HIGGS.csv1.csv", header = F)
    }
    else{
      df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1.csv", header = F)
    }
    
    #samplesize = dim(df)[1]
    samplesize = 10000
    X = cbind(rep(1, samplesize), df[1:samplesize,23:29])
    colnames(X)[1] <- "V1"
    y = df[1:samplesize,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs1_large"){
    # subset of the higgs data set, 10**5 lines
    if(server){
      df <- read.csv("/mrc-bsu/scratch/ab2603/higgsdata/HIGGS.csv1_large.csv", header = F)
    }
    else{
      df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1_large.csv", header = F)
    }
    
    samplesize = dim(df)[1]
    X = cbind(rep(1, samplesize), df[,23:29])
    colnames(X)[1] <- "V1"
    y = df[,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs2_large"){
    # subset of the higgs data set, 10**5 lines
    if(server){
      df <- read.csv("/mrc-bsu/scratch/ab2603/higgsdata/HIGGS.csv1_large.csv", header = F)
    }
    else{
      df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1_large.csv", header = F)
    }
    
    samplesize = dim(df)[1]
    X = cbind(rep(1, samplesize), df[,2:22])
    colnames(X)[1] <- "V1"
    y = df[,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs2_small"){
    # subset of the higgs data set, 10**5 lines
    if(server){
      df <- read.csv("/mrc-bsu/scratch/ab2603/higgsdata/HIGGS.csv1.csv", header = F)
    }
    else{
      df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1.csv", header = F)
    }
    
    samplesize = 10000
    X = cbind(rep(1, samplesize), df[,2:22])
    colnames(X)[1] <- "V1"
    y = df[,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs2_full"){
    # subset of the higgs data set
    if(server){
      df <- read.csv("/mrc-bsu/scratch/ab2603/higgsdata/HIGGS.csv", header = F)
    }
    else{
      df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv", header = F)
    }
    
    samplesize = dim(df)[1]
    X = cbind(rep(1, samplesize), df[,2:22])
    colnames(X)[1] <- "V1"
    y = df[,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs2"){
    # subset of the higgs data set
    df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1.csv", header = F)
    samplesize = dim(df)[1]
    X = cbind(rep(1, samplesize), df[,2:22])
    colnames(X)[1] <- "V1"
    y = df[,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs3"){
    # subset of the higgs data set
    df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1.csv", header = F)
    samplesize = dim(df)[1]
    X = cbind(rep(1, samplesize), df[,2:29])
    colnames(X)[1] <- "V1"
    y = df[,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "sim1"){
    # this is a missspecified model, the last observation is missing
    set.seed(123)
    samplesize = nobs
    betasize = 9
    
    X = matrix(rnorm(n = samplesize*betasize), ncol = betasize)
    colnames(X) <- paste("X", c(1:betasize), sep="")
    betacoef = c(-1,1,rep(0, betasize-3),1)
    prodbeta = X%*%betacoef
    y = rbinom(samplesize, size= 1, prob=1/(1+exp(-prodbeta)))
    list_data[["X"]] <- X[,1:(betasize-1)] # remove the last observation here
    list_data[["y"]] <- y
    list_data[["betastar"]] <- betacoef[1:(betasize-1)]
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "sim2"){
    set.seed(123)
    samplesize = nobs
    betasize = 5
    if(highcorr) {
      corr_mat <- matrix(0.99, nrow = betasize, ncol = betasize)
      diag(corr_mat) <- rep(1,betasize)
      
      var_x <- rep(1,betasize)
      sigma_x <- diag(var_x) %*% corr_mat %*% diag(var_x)
    }
    else sigma_x <- diag(nrow = betasize)
    X = matrix(rnorm(n = samplesize*betasize), ncol = betasize)
    for (n in 1:samplesize) {
      X[n,] <- mvtnorm::rmvnorm(1, sigma=sigma_x)
    }
    colnames(X) <- paste("X", c(1:betasize), sep="")
    betacoef1 = c(-1,1,rep(0, betasize-3),1)
    betacoef2 = c(-1,1,0.01, rep(0, betasize-4),1)
    #browser()
    prodbeta1 = X[1:as.integer(samplesize/2),]%*%betacoef1
    prodbeta2 = X[(as.integer(samplesize/2)+1):samplesize,]%*%betacoef2
    prodbeta = c(prodbeta1, prodbeta2)
    y = rbinom(samplesize, size= 1, prob=1/(1+exp(-prodbeta)))
    list_data[["X"]] <- X
    list_data[["y"]] <- y
    list_data[["betastar"]] <- list(betacoef1, betacoef2)
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "flights_simple"){
    df <- flights %>%  mutate(y = arr_delay>0) %>% mutate(intercept = 1)%>% select(y, intercept, dep_delay) %>% na.omit()
    df <- data.matrix(df)
    X = df[,2:3]
    y = df[,1]*1
    # lm <- glm(y ~  dep_delay, data = df, family = "binomial")
    # summary(lm)
    list_data[["X"]] <- X
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "flights_complex1"){
    df <- flights %>%  mutate(y = arr_delay>0) %>% select(y, dep_delay, carrier) %>% na.omit()
    df <- df %>% to_dummy(carrier, suffix = "label") %>%
      bind_cols(df) %>%
      select(y, dep_delay, everything())
    #carrier_list <- unique(df$carrier)
    #for(carrier_i in carrier_list){
    #  df <- df %>% mutate()
    #}
    df <- data.matrix(df)
    X = df[,2:18]
    y = df[,1]*1
    list_data[["X"]] <- X
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "flights_complex2"){
    print("run the complex model with independent intercepts")
    df <- flights %>%  mutate(y = arr_delay>0) %>% select(y, dep_delay, carrier) %>% na.omit()
    df <- df %>% to_dummy(carrier, suffix = "label") %>%
      bind_cols(df) %>%
      select(y, dep_delay, everything())
    
    carrier_list <- unique(df$carrier)
    for(carrier_i in carrier_list){
      df[[paste("delay_", carrier_i, sep="")]] <- df[[paste("carrier_", carrier_i, sep="")]]*df[["dep_delay"]]
    }
    df %<>% select(-carrier)
    df %<>% select(-dep_delay)
    num_covars <- dim(df)[2]
    df <- data.matrix(df)
    X = df[,2:num_covars]
    y = df[,1]*1
    list_data[["X"]] <- X
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else{print("The type of dataset does not exist!")}
  return(list_data)
}
