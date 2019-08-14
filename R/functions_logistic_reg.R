# functions logistic regression
#library(mvtnorm)
#library(boot)
#library(bridgesampling)
#library(MASS)
#library(car)
#library(BayesLogit)
#source("R/logisticreg/functions_product_normal.R")

dMvn <- function(X,mean,Sigma,log=F) {
  k <- ncol(X)
  rooti <- backsolve(chol(Sigma),diag(k))
  quads <- colSums((crossprod(rooti,(t(X)-mean)))^2)
  logres <- -(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads
  if(log==T){
    return(exp(logres))
  }
  else{
    return(logres)
  }

}

f_rnorm <- function(mu, Sigma){
  d <- length(mu)
  risotropic <- rnorm(d)
  uppertri <- chol(Sigma)
  res <- t(uppertri) %*% risotropic + mu
  return(res)
}

f_dmvnorm <-  function(x, mu, Sigma, logflag = F){
  # own implementation of the dmvnorm
  Sigma_inv <- solve(Sigma)
  matprod <- -0.5*t(x - mu) %*% Sigma_inv %*% (x - mu)
  lognormconst <- -0.5*log(det(2*pi*Sigma))
  if(logflag == F){
    res <- exp(matprod + lognormconst)
  }
  else{
    res <- matprod + lognormconst
  }
  return(res)
}


f_prep_prior_logistic <- function(ssplitlist, scale=1){
  # uses a normal prior
  #browser()
  ssplits = length(ssplitlist)
  d = ssplitlist[[1]][["Pparams"]]
  mu = rep(0, d)
  Sigma = ssplits*diag(rep(scale,d), nrow=d)
  for(isplit in 1:ssplits){
    ssplitlist[[isplit]][["bprior"]] <- mu
    ssplitlist[[isplit]][["scale"]] <- scale*ssplits
    ssplitlist[[isplit]][["Bprior"]] <- Sigma
    ssplitlist[[isplit]][["Bpriorinv"]] <- solve(Sigma)
    ssplitlist[[isplit]][["ssplits"]] <- ssplits
    ssplitlist[[isplit]][["isplit"]] <- isplit
  }
  return(ssplitlist)
}


# functions for logistic regression
f_pg_cond <- function(betacurrent, params){
  # a function that samples from the polya gamma sampler
  #browser()
  parapg <- params[["X"]] %*% betacurrent
  #browser()
  #z <- apply(parapg, 1, function(x) pgdraw(1,x))
  #z <- BayesLogit::rpg(num = params$n, h=1, z = parapg)
  z <- pgdraw(1, parapg)
  return(z)
}

f_norm_cond <- function(Zsamples, params){
  # a function that samples from the conditional normal
  #browser()
  # old version
  #Omega <- diag(Zsamples)
  #Vw <- solve(t(params[["X"]])%*%Omega%*%params[["X"]]+params[["Bpriorinv"]])
  # new version
  intermatprod <-t(params[["X"]]*Zsamples)%*%params[["X"]]+params[["Bpriorinv"]]
  Vw <- solve(intermatprod)
  Vw <- (t(Vw)+Vw)/2
  mw <- Vw%*%(t(params[["X"]])%*%params[["kappa"]]+params[["Bpriorinv"]]%*%params[["bprior"]])
  #beta <- mvrnorm(1, mu = mw, Sigma = Vw)
  beta <- f_rnorm(mu = mw, Sigma = Vw)
  return(list(beta=beta, mw=mw, Vw=Vw))
}

f_sample_chain <- function(n_steps, params, burnin=500){
  # bring the gibbs sampler to run
  # prepare the data
  res_list <- list()
  res_list[["chainbeta"]] <- matrix(0, nrow=params[["d"]], ncol = n_steps)
  #res_list[["chainz"]] <- matrix(0, nrow=params[["n"]], ncol = n_steps)
  res_list[["chainmw"]] <- matrix(0, nrow=params[["d"]], ncol = n_steps)
  res_list[["chainVw"]] <- array(0, dim=c(params[["d"]], params[["d"]], n_steps))
  #res_list[["chainLambda"]] <- array(0, dim=c(params[["d"]], params[["d"]], n_steps))
  #res_list[["chaineta"]] <- array(0, dim=c(params[["d"]], n_steps))
  #res_list[["chainxi"]] <- array(0, dim=c(1, n_steps))

  # samples first step
  betacurrent = mvrnorm(n = 1, mu = params[["bprior"]], Sigma = params[["Bprior"]])
  res_list[["chainbeta"]][,1] <- betacurrent
  #res_list[["chainz"]][,1] <- f_pg_cond(res_list[["chainbeta"]][,1], params)
  #chainz <- f_pg_cond(res_list[["chainbeta"]][,1], params)
  res_list[["chainmw"]][,1] <- params[["bprior"]]
  res_list[["chainVw"]][,,1] <- params[["Bprior"]]
  #res_list[["chainLambda"]][,,1] <- solve(params[["Bprior"]])
  #browser()
  #res_list[["chaineta"]][,1] <- res_list[["chainLambda"]][,,1]%*%res_list[["chainmw"]][,1]
  #res_list[["chainxi"]][,1] <- -0.5*(params$d*log(2*pi) - log(det(res_list[["chainLambda"]][,,1])) +
  #                                     t(res_list[["chaineta"]][,1])%*%res_list[["chainLambda"]][,,1]%*%res_list[["chaineta"]][,1])
  # loop over the iterations
  for(tsteps in 2:n_steps){
    #res_list[["chainz"]][,tsteps] <- f_pg_cond(res_list[["chainbeta"]][,tsteps-1], params)
    #browser()
    chainz <- f_pg_cond(res_list[["chainbeta"]][,tsteps-1], params)
    interres <- f_norm_cond(chainz, params)
    res_list[["chainbeta"]][,tsteps] <- interres$beta
    res_list[["chainmw"]][,tsteps] <- interres$mw
    res_list[["chainVw"]][,,tsteps] <- interres$Vw
    #res_list[["chainLambda"]][,,tsteps] <- solve(interres$Vw)
    #browser()
    #res_list[["chaineta"]][,tsteps] <- res_list[["chainLambda"]][,,tsteps]%*%res_list[["chainmw"]][,tsteps]
    #res_list[["chainxi"]][,tsteps] <- f_help_closed_form_normal(params$d, res_list[["chainLambda"]][,,tsteps], res_list[["chaineta"]][,tsteps]) #-0.5*(params$d*log(2*pi) - log(det(res_list[["chainLambda"]][,,tsteps])) +
    #t(res_list[["chaineta"]][,tsteps])%*%res_list[["chainLambda"]][,,tsteps]%*%res_list[["chaineta"]][,tsteps])
  }
  res_list[["chainbeta"]] <- res_list[["chainbeta"]][,burnin:n_steps]
  #res_list[["chainz"]] <- res_list[["chainz"]][,burnin:n_steps]
  res_list[["chainmw"]] <- res_list[["chainmw"]][,burnin:n_steps]
  res_list[["chainVw"]] <- res_list[["chainVw"]][,,burnin:n_steps]
  #res_list[["chainLambda"]] <- res_list[["chainLambda"]][,,burnin:n_steps]
  #res_list[["chaineta"]] <- res_list[["chaineta"]][,burnin:n_steps]
  #res_list[["chainxi"]] <- res_list[["chainxi"]][,burnin:n_steps]
  return(res_list)
}

f_help_closed_form_normal <- function(d, Lambda, eta){
  res <- -0.5*(d*log(2*pi) - log(det(Lambda)) + t(eta)%*%Lambda%*%eta)
  return(res)
}

f_loglik <- function(betastar, params){
  # log likelihood of the logistic regression
  if(!is.null(dim(betastar))){
    il_xbetaprod = inv.logit(params[["X"]]%*%t(betastar))
  }
  else{
    il_xbetaprod = inv.logit(params[["X"]]%*%betastar)
  }
  part1 = il_xbetaprod**params[["y"]]
  part2 = (1-il_xbetaprod)**(1-params[["y"]])
  return(apply(log(part1)+log(part2), 2, sum))
}

f_logprior <- function(x, params){
  # log likelihood of prior
  return(dmvnorm(x, mean= params[["bprior"]], sigma=params[["Bprior"]], log = TRUE))
  #browser()
  #return(dmvnrm_arma(matrix(x, ncol = params$d), mean= params[["bprior"]], sigma=params[["Bprior"]], logd = TRUE))
  #return(dMvn(matrix(x, ncol = params$d), mean = params[["bprior"]], Sigma=params[["Bprior"]], log = TRUE))

}

f_rprior <- function(n, params){
  # log likelihood of prior
  return(mvrnorm(n, mu= params[["bprior"]], Sigma=params[["Bprior"]]))
}


f_logconditional_mean <- function(x, reschain, params){
  # log likelihood of prior
  inter = 0
  for(i in c(1:dim(reschain$chainbeta)[2])){
    beta = reschain$chainmw[,i]
    sigma = reschain$chainVw[,,i]
    inter = inter+dmvnorm(x, mean=beta, sigma=sigma, log=FALSE)
    #browser()
    #inter = inter+exp(dmvnrm_arma(matrix(x, ncol=params$d), mean=beta, sigma=sigma, logd=T))
    #inter = inter+dMvn(matrix(x, ncol=params$d), mean=beta, Sigma=sigma, log=F)
  }
  return(log(inter/i))
}

f_logposterior <- function(x, params.=params, data=NULL) f_loglik(x, params.)+f_logprior(x, params.)

log_sum_exp <- function(x){
  # numerically stable version of log(sum(exp(x)))
  log(sum(exp(x - max(x)))) + max(x)
}

f_parallel_sampler <- function(datasplits, n_steps=22000, burnin=1000, filename_individual_sim = "tempsavesim", returnres=F){
  # run the sampler on the datasplits
  ssplits <- length(datasplits)
  #browser()
  list_rep_chain <- list()
  for(i_s in c(1:ssplits)){
    # run the sampler here
    print(paste("sample split ", i_s))
    srepchain <- f_sample_chain(n_steps=n_steps, params=datasplits[[i_s]], burnin=burnin)
    #browser()
    muproxy <- apply(srepchain$chainbeta, 1, mean)
    normconstchib <- f_logprior(muproxy, datasplits[[i_s]])-
      f_logconditional_mean(muproxy, srepchain, datasplits[[i_s]])+
      f_loglik(muproxy, datasplits[[i_s]])
    srepchain[["chib"]] <- normconstchib
    print(normconstchib)
    # lb <- rep(-Inf, datasplits[[i_s]]$d)
    # ub <- rep(Inf, datasplits[[i_s]]$d)
    # betaposterior <- data.frame(t(srepchain$chainbeta))
    # colnames(betaposterior) <- colnames(datasplits[[i_s]][["X"]])
    # names(lb) <- names(ub) <- colnames(betaposterior)
    # bridge_result <- bridge_sampler(samples = as.matrix(betaposterior),
    #                                 data = NULL,
    #                                 log_posterior = f_logposterior,
    #                                 params =  datasplits[[i_s]],
    #                                 lb = lb, ub = ub, silent = TRUE)
    # srepchain[["bridge"]] <- bridge_result$logml
    save(srepchain, file = paste(filename_individual_sim, toString(i_s), ".RData", sep = ""))
    if(returnres==F) rm(srepchain)
    else list_rep_chain[[i_s]] <- srepchain
  }
  if(returnres==T) return(list_rep_chain)
  else return(list())
}

# functions daniel
normalise_product_norm <- function(mu, Sigma){
  # return the log of normalising constant for a product of s normal distributions
  # each normal distribution is d-dimensional
  # mu is a d times s matrix where each column gives the mean vector
  # Sigma is an array of dimension (d,d,s) where each slice (,,i) gives a covariance matrix
  s <- ncol(mu)
  d <- nrow(mu)
  Lambda <- array(0, dim=dim(Sigma))
  eta <- array(0, dim=dim(mu))
  for (i in 1:s){
    Lambda[,,i] <- solve(Sigma[,,i])
    eta[,i] <- Lambda[,,i] %*% mu[,i]
  }
  etasum <- rowSums(eta)
  Lambdasum <- as.matrix(apply(Lambda, c(1,2), sum))
  sumdet <- sum(vapply(1:s, function(i) log(det(as.matrix(Lambda[,,i]))), matrix(0,1,1)))
  sumquad <- sum(vapply(1:s, function(i) t(eta[,i]) %*% Sigma[,,i] %*% eta[,i], matrix(0,1,1)))
  Z <- (-1/2*s*d*log(2*base::pi)+1/2*sumdet-1/2*sumquad)-(-1/2*d*log(2*base::pi)+1/2*log(det(Lambdasum))-1/2*t(etasum)%*%solve(Lambdasum)%*%etasum)
  return(Z)
}

# suprior constant
if(F){
  d <- 5
  ssplits <- 3
  n <- 10**1
  Vprior <- diag(rep(1, d))
  sVprior <- ssplits*Vprior
  is_samples <- rmvnorm(n, mean = rep(0,d), sigma = sVprior)
  logweights <- (1/ssplits)*dmvnorm(is_samples, sigma = Vprior, log=T)-dmvnorm(is_samples, sigma = sVprior, log=T)
}

f_alpha_sub <- function(ssplits, Vprior, typeprior="normal"){
  if(typeprior=="normal"){
    alpha1 <- det((2*pi)*Vprior)**(-1/(2*ssplits))
    alpha2 <- det((2*pi)*ssplits*Vprior)**(1/2)
    alphasubpriorconstant <- (log(alpha1) + log(alpha2))
  }
  else if (typeprior=="laplace_normal"){
    # assuming only the last two components are normal
    d <- dim(Vprior)[1]
    # the normal part
    #browser()
    alpha1 <- det((2*pi)*Vprior[(d-1):d, (d-1):d])**(-1/(2*ssplits))
    alpha2 <- det((2*pi)*ssplits*Vprior[(d-1):d, (d-1):d])**(1/2)
    part1 <- (log(alpha1)+log(alpha2))
    # the expoential part, assuming that all the prior variance are the same
    part2 <- (-((d - 2)/ssplits)*log(2*Vprior[1,1])) + ((d - 2)*log(2*ssplits*Vprior[1,1]))
    alphasubpriorconstant <- part1+part2
  }
  else if (typeprior=="laplace"){
    # assuming everything is a laplace prior
    d <- dim(Vprior)[1]
    # the expoential part, assuming that all the prior variance are the same
    part2 <- (-((d )/ssplits)*log(2*Vprior[1,1])) + ((d )*log(2*ssplits*Vprior[1,1]))
    alphasubpriorconstant <- part2
  }
  return(alphasubpriorconstant)
}

# this is our own function
logistic_isub_adapted <- function(filename_individual_sim, ssplits){
  # models is a list of length s where each element is the augmented output of a BayesLogit fit.
  # augmented logitfit includes conditional mean and conditional variance of beta at each iteration.
  # each Gibbs sampler run needs to be the same number of iterations (B)
  load(paste(filename_individual_sim, toString(1), ".RData", sep = ""))
  #browser()
  p <- ncol(t(srepchain$chainbeta))
  B <- nrow(t(srepchain$chainbeta))
  condmu <- array(0, dim = c(ssplits,p,B))
  condSigma <- array(0, dim = c(ssplits, p,p,B))
  for (i in 1:ssplits){
    load(paste(filename_individual_sim, toString(i), ".RData", sep = ""))
    condmu[i,,] <- srepchain$chainmw
    condSigma[i,,,] <- srepchain$chainVw
  }
  #browser()
  if(ssplits == 1){
    print("code does not work for a single split")
  }
  #f_integral_product_gaussian(condmu[,,1], aperm(condSigma[,,,1], c(2,3,1)))
  condintegral <- vapply(1:B, function (b) f_integral_product_gaussian(condmu[,,b], aperm(condSigma[,,,b], c(2,3,1))), matrix(0,1,1))
  #condintegral <- sapply(1:B, function (b) normalise_product_norm(t(condmu[,,b]),aperm(condSigma[,,,b], c(2,3,1))))
  logIsub <- log_sum_exp(condintegral) - log(B)
  ESS <- Ess_logweights(condintegral)*100/B
  print(paste("ESS", ESS))
  reslist <- list()
  reslist[["logIsub"]] <- logIsub
  reslist[["ESS"]] <- ESS
  reslist[["condintegral"]] <- condintegral
  return(reslist)
}

Ess_logweights <- function(logweights){
  # calculates the ESS based on the log weights
  numerator <- 2*log_sum_exp(logweights)
  denominator <- log_sum_exp(2*logweights)
  return(exp(numerator - denominator))
}

augment_logitfit <- function(logitfit, m0, V0){
  # take output of BayesLogit function and calculate parameters of conditional normal distribution at each iteration
  # m0 is prior mean
  # V0 is prior variance. BayesLogit takes prior precision matrix.

  P0 <- solve(V0)
  B <- nrow(logitfit$beta)
  y <- logitfit$y
  X <- logitfit$X
  n <- logitfit$n
  p <- ncol(logitfit$beta)
  beta_mu <- array(0, dim=c(p, B))
  beta_Sigma <- array(0, dim=c(p, p, B))
  kappa <- y*n-n/2
  #browser()
  for (b in 1:B){
    w <- as.vector(logitfit$w[b,])
    Vw <- solve(car::wcrossprod(X, X, w)+P0)
    Vw <- (t(Vw)+Vw)/2
    mw <- Vw %*% (t(X) %*% kappa+(P0)%*%m0)
    beta_mu[,b] <- mw
    beta_Sigma[,,b] <- Vw
    # beta_mu[,b] gives the conditional mean of beta at iteration b in the Gibbs run
    # beta_Sigma[,,b] gives the conditional variance of beta at iteration b in the Gibbs run
  }
  logitfit$beta_mu <- beta_mu
  logitfit$beta_Sigma <- beta_Sigma
  return(logitfit)
}

chib_logitfit <- function(logitfit, m0, V0, betastar){
  # use Chib's method to estimate log evidence given logitfit object
  # betastar is test point to plug into Candidate's formula
  # m0 is prior mean and V0 is prior variance
  if (is.null(logitfit$beta_mu)){
    logitfit <- augment_logitfit(logitfit, m0, V0)
  }
  P0 <- solve(V0)
  betasamp <- logitfit$beta
  B <- nrow(betasamp)
  b <- ncol(betasamp)
  betastarmat <- matrix(betastar, ncol=b)
  y <- logitfit$y
  X <- logitfit$X
  n <- logitfit$n
  m <- logitfit$beta_mu
  V <- logitfit$beta_Sigma
  logdens <- numeric(B)
  for (b in 1:B){
    logdens[b] <- mvtnorm::dmvnorm(betastar, mean=m[,b], sigma=V[,,b], log=TRUE)
    #logdens[b] <- dmvnrm_arma(betastarmat, mean=m[,b], sigma=V[,,b], logd=TRUE)
    #logdens[b] <- dMvn(betastarmat, mean=m[,b], Sigma=V[,,b], log=TRUE)
  }
  loglk <- lrlk(y, X, betastar, n)
  logprior <-  mvtnorm::dmvnorm(betastar,mean=m0, sigma=V0, log=TRUE)
  logpost <- log_sum_exp(logdens)-log(B)
  evidence <- loglk+logprior-logpost
  return(evidence)
}

log_sum_exp <- function(x){
  # numerically stable version of log(sum(exp(x)))
  log(sum(exp(x - max(x)))) + max(x)
}

log_mean_exp <- function(x){
  # numerically stable version of log(sum(exp(x)))
  log(mean(exp(x - max(x)))) + max(x)
}


lrlk <- function(y, X, beta, n=rep(1, length(y))){
  # logistic regression log likelihood
  # y is vector of length N
  # X is N times p design matrix
  # beta is p-length vector
  # y is vector of outcomes (normally 0,1) but can be success proportions if there are duplicate covariates in design matrix
  # see BayesLogit logit.combine function
  # n is vector of counts, necessary if using summary stat representation of full dataset
  eta <- X %*% beta
  sum(beta*(t(X) %*% (y*n)))-sum(n*log(1+exp(eta)))
}

logistic_isub <- function(models){
  # models is a list of length s where each element is the augmented output of a BayesLogit fit.
  # augmented logitfit includes conditional mean and conditional variance of beta at each iteration.
  # each Gibbs sampler run needs to be the same number of iterations (B)
  s <- length(models)
  p <- ncol(models[[1]]$beta)
  B <- nrow(models[[1]]$beta)
  condmu <- array(0, dim=c(s,p,B))
  condSigma <- array(0, dim=c(s, p,p,B))
  for (i in 1:s){
    condmu[i,,] <- models[[i]]$beta_mu
    condSigma[i,,,] <- models[[i]]$beta_Sigma
  }
  #browser()
  condintegral <- vapply(1:B, function (b) normalise_product_norm(t(condmu[,,b]),aperm(condSigma[,,,b], c(2,3,1))), matrix(0,1,1))
  logIsub <- log_sum_exp(condintegral)-log(B)
  ESS <- Ess_logweights(condintegral)*100/B
  print(paste("ESS", ESS))
  reslist <- list()
  reslist[["logIsub"]] <- logIsub
  reslist[["ESS"]] <- ESS
  return(reslist)
}
