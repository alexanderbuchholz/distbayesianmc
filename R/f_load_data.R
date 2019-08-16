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
f_dataset_loader <- function(dataset="pima", nobs=5*10**3, highcorr = T){
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
    df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv", header = F)
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
    df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1.csv", header = F)
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
    df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1_large.csv", header = F)
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
    df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv1_large.csv", header = F)
    samplesize = dim(df)[1]
    X = cbind(rep(1, samplesize), df[,2:22])
    colnames(X)[1] <- "V1"
    y = df[,1]
    list_data[["X"]] <- X # remove the last observation here
    list_data[["y"]] <- y
    list_data[["dataset"]] <- dataset
  }
  else if (dataset == "higgs2_full"){
    # subset of the higgs data set
    df <- read.csv("/scratch/alexander/higgsdata/HIGGS.csv", header = F)
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
