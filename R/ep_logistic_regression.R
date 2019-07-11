# logistic regression using ep
#library(EPGLM)
#library(BayesLogit)
#library(MASS)

if(FALSE){
  data_pima = Pima.tr
  X = scale(data_pima[,1:7])
  ones = rep(1, nrow(data_pima))
  X = cbind(ones, X)
  y = (data_pima[,8] == "Yes")*1
  Sol<-EPlogit(X,y,1)

  source("R/helperfunctions/f_helper.R")
  source("R/logisticreg/functions_logistic_reg.R")
  ssplits = 4
  datasplits <- f_pack_split_data(X,y, ssplits = ssplits)
  scale = 1
  datasplits <- f_prep_prior_logistic(datasplits, scale=scale)
}

#Sol<-EPlogit(X, y, scale)


#covarproxy <- Sol$V
#muproxy <- Sol$m
#samples_is <- rmvnorm(n= 20000, mean= muproxy, sigma=covarproxy)
#logweights <- f_loglik(samples_is, params = datasplits[[1]])+f_logprior(samples_is, params = datasplits[[1]]) - dmvnorm(samples_is, mean=muproxy, sigma=covarproxy, log=TRUE)
#normconstis <- log(mean(exp(logweights)))
#Sol$Z

f_ep_splits <- function(datasplits){
  # apply ep
  eplist <- list()
  icounter <- 1
  r_joint <- 0
  Q_joint <- 0
  n <- 10**5
  # loop over data set
  for(datasp in datasplits){
    Sol<-EPlogit(datasp$X, datasp$y, datasp$scale)
    eplist[[icounter]] <- list()
    eplist[[icounter]][["m_ep"]] <- Sol$m
    eplist[[icounter]][["V_ep"]] <- Sol$V
    eplist[[icounter]][["Z_ep"]] <- Sol$Z
    eplist[[icounter]][["Q_ep"]] <- solve(eplist[[icounter]][["V_ep"]])
    eplist[[icounter]][["r_ep"]] <- eplist[[icounter]][["Q_ep"]]%*%Sol$m
    r_joint <- r_joint+eplist[[icounter]][["r_ep"]]
    Q_joint <- Q_joint+eplist[[icounter]][["Q_ep"]]

    # calculate exact log Z using IS
    samples_is <- mvrnorm(n, mu = Sol$m, Sigma = Sol$V)
    logweightsinter <- f_loglik(samples_is, params = datasp)+f_logprior(samples_is, params = datasp) - dmvnorm(samples_is, mean=Sol$m, sigma=Sol$V, log=TRUE)
    eplist[[icounter]][["Z_exact"]] <- log_sum_exp(logweightsinter)-log(n)

    icounter <- icounter+1
  }
  eplist[["r_joint"]] <- r_joint
  eplist[["Q_joint"]] <- Q_joint
  eplist[["Sigma_joint"]] <- solve(Q_joint)
  eplist[["mu_joint"]] <- eplist[["Sigma_joint"]]%*%r_joint
  return(eplist)
}



f_logistic_splits <- function(datasplits){
  # apply ep
  logisticlist <- list()
  icounter <- 1
  r_joint <- 0
  Q_joint <- 0
  n <- 10**5
  # loop over data set
  for(datasp in datasplits){
    data_inter <- as.data.frame(cbind(datasp$y, datasp$X))
    pvars <- dim(data_inter)[2]
    dfnames <- colnames(data_inter)
    colnames(data_inter)[2:pvars] <- paste("X", dfnames[2:pvars], sep="")

    form <- as.formula(paste(dfnames[1], "~", paste(c("0",colnames(data_inter)[2:pvars]), collapse="+")))
    res <- glm(form, data=data_inter, family = "binomial")
    res$coefficients
    logisticlist[[icounter]] <- list()
    logisticlist[[icounter]] <- res$coefficients
    # calculate exact log Z using IS
    icounter <- icounter+1
  }
  return(logisticlist)
}




f_ep_splitted_is <- function(eplist, datasplits, n=10**4){
  # this function is not useful
  ssplits = length(datasplits)
  samples_is <- mvrnorm(n, mu = eplist$mu_joint, Sigma = eplist$Sigma_joint)
  logweights <- 0
  internormconstant <- 0
  for(i in 1:ssplits){
    logweightsinter <- f_loglik(samples_is, params = datasplits[[i]])+f_logprior(samples_is, params = datasplits[[i]]) - dmvnorm(samples_is, mean=eplist[[i]][["m_ep"]], sigma=eplist[[i]][["V_ep"]], log=TRUE)
    logweights <- logweights+logweightsinter
    internormconstant <-internormconstant + log_sum_exp(logweightsinter)-log(n)
  }
  ESS = (sum(exp(logweights))**2)/sum(exp(logweights)**2)
  print(ESS*100/n)
  return(log_sum_exp(logweights)-log(n)-internormconstant)
}

#eplist <-f_ep_splits(datasplits)
#boxplot(sapply(1:100, function(i){f_ep_splitted_is(eplist, datasplits, n=10**3)}))

