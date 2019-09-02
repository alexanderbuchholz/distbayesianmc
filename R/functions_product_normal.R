# is product normal
#library(mvtnorm)
#library(stats)
#source("R/logisticreg/functions_logistic_reg.R")
if(F){
  set.seed(42)

  d <- 2
  ssplits <- 3
  mu_all <- rmvnorm(ssplits, mean = rep(0, d), sigma = diag(rep(1, d)))
  Sigma_all <- rWishart(n = ssplits, df = 5, Sigma = diag(rep(1, d)))

  n <- 10**6
  is_samples <- rmvnorm(n, mean = mu_all[1,], sigma = Sigma_all[,,1])
  log_weights <- dmvnorm(is_samples, mean = mu_all[2,], sigma = Sigma_all[,,2], log = T)+dmvnorm(is_samples, mean = mu_all[3,], sigma = Sigma_all[,,3], log = T)
  is_approx1 = log_sum_exp(log_weights) - log(n)
  is_approx1
  #
  is_samples <- rmvnorm(n, mean = mu_all[2,], sigma = Sigma_all[,,2])
  log_weights <- dmvnorm(is_samples, mean = mu_all[1,], sigma = Sigma_all[,,1], log = T)+dmvnorm(is_samples, mean = mu_all[3,], sigma = Sigma_all[,,3], log = T)
  is_approx2 = log_sum_exp(log_weights) - log(n)
  is_approx2

  is_samples <- rmvnorm(n, mean = mu_all[3,], sigma = Sigma_all[,,3])
  log_weights <- dmvnorm(is_samples, mean = mu_all[1,], sigma = Sigma_all[,,2], log = T)+dmvnorm(is_samples, mean = mu_all[3,], sigma = Sigma_all[,,3], log = T)
  is_approx3 = log_sum_exp(log_weights) - log(n)
  is_approx3

}

f_gaussian_standard_to_natural_params <- function(mu_all, Sigma_all){
  # transforms standard gaussian parameters to natural parameters
  ssplits <- dim(Sigma_all)[3]
  d <- dim(Sigma_all)[1]
  #browser()
  Lambda_all <-tryCatch({
    #Lambda_all <- 
    Lambda_all <- array(unlist(lapply(c(1:ssplits), function(i) solve(Sigma_all[,,i]))), dim=c(d,d,ssplits))
    eta_all <- t(sapply(c(1:ssplits), function(i) Lambda_all[,,i]%*%mu_all[i,]))
    res_list <- list()
    res_list[["eta_all"]] <- eta_all
    res_list[["Lambda_all"]] <- Lambda_all
    return(res_list)
  }, error=function(...){
      #browser()
      Lambda_all <- array(unlist(lapply(c(1:ssplits), function(i) solve(Sigma_all[,,i]+diag(x=10^(-6), d, d)))), dim=c(d,d,ssplits))
      eta_all <- t(sapply(c(1:ssplits), function(i) Lambda_all[,,i]%*%mu_all[i,]))
      res_list <- list()
      res_list[["eta_all"]] <- eta_all
      res_list[["Lambda_all"]] <- Lambda_all
      return(res_list)
  })
  
  
  
}

if(F){
  res_list <- f_gaussian_standard_to_natural_params(mu_all, Sigma_all)
  res_list_inverse <- f_gaussian_standard_to_natural_params(res_list$eta_all, res_list$Lambda_all)

  sum((res_list_inverse$eta_all-mu_all)**2)
  sum((res_list_inverse$Lambda_all-Sigma_all)*2)

}

f_log_norm_const_natural_gaussian_standard_params <- function(mu, Sigma){
  # the parameters are with respect to a standard gaussian,
  # not the natural parameters
  Sigma_inv <- solve(Sigma)
  res <- 0.5*mu%*%Sigma_inv%*%mu+0.5*log(det(2*pi*Sigma))
  return(res)
}

f_log_norm_const_natural_gaussian_natural_params <- function(eta, Lambda){
  # the parameters are with respect to the natural parameters
  #browser()
  Lambda_inv <- solve(Lambda)
  #logdet <- log(det(2*pi*Lambda_inv))
  logdet <- unlist(determinant(2*pi*Lambda_inv, logarithm = T))[1]
  #res <- 0.5*eta%*%Lambda_inv%*%eta+0.5*log(det(2*pi*Lambda_inv))
  res <- 0.5*eta%*%Lambda_inv%*%eta+0.5*logdet
  return(res)
}

if(F){
  f_log_norm_const_natural_gaussian_standard_params(mu_all[1,], Sigma_all[,,1])
  f_log_norm_const_natural_gaussian_natural_params(res_list$eta_all[1,], res_list$Lambda_all[,,1])
}

f_log_norm_const_gaussian <- function(x, mu, Sigma){
  Sigma_inv <- solve(Sigma)
  part1 <- -0.5*t(x-mu)%*%Sigma_inv%*%(x-mu)
  part2 <- -0.5*log(det(2*pi*Sigma))
  return(part1+part2)
}

f_exponential_form_normal <- function(x, mu, Sigma){
  # without normalizing constant
  Sigma_inv <- solve(Sigma)
  #browser()
  res <- (-0.5*(x%*%Sigma_inv%*%x))+(mu%*%Sigma_inv%*%x)
  return(res)
}
if(F){
  n = 10**1
  normal_is <- rmvnorm(n, mean=mu_all[1,], sigma = Sigma_all[,,1])

  logweights <- sapply(c(1:n), function(i) {f_exponential_form_normal(normal_is[i,], mu_all[1,], Sigma_all[,,1])-dmvnorm(normal_is[i,], mean=mu_all[1,], sigma = Sigma_all[,,1], log = T)})
  log_sum_exp(logweights) -log(n)


  f_exponential_form_normal(normal_is[1,], mu_all[1,], Sigma_all[,,1])-
    f_log_norm_const_natural_gaussian_standard_params(mu_all[1,], Sigma_all[,,1])-
    dmvnorm(normal_is[1,], mean=mu_all[1,], sigma = Sigma_all[,,1], log = T)
}
#logweights <- sapply(c(1:n), function(i) f_exponential_form_normal(normal_is[i,], mu_all[1,], diag(rep(1,d)))-dmvnorm(normal_is[i,], mean=mu_all[1,], sigma = diag(rep(1,d)), log = T))

f_sum_params_gaussian <- function(mu_all, Sigma_all){
  ssplits <- dim(mu_all)[1]
  res_list <- f_gaussian_standard_to_natural_params(mu_all, Sigma_all)
  #Q_joint <- rowSums(res_list$Lambda_all, dims = 2)
  Q_joint <- apply(res_list$Lambda_all, c(1,2), sum)
  r_joint <- apply(res_list$eta_all, 2, sum)
  Sigma_joint <- solve(Q_joint)
  mu_joint <- Sigma_joint %*% r_joint
  res <- list()
  res[["mu_joint"]] <-  mu_joint
  res[["Sigma_joint"]] <-  Sigma_joint
  return(res)
}

f_integral_product_gaussian <- function(mu_all, Sigma_all){
  ssplits <- dim(mu_all)[1]
  #browser()
  res_list <- f_gaussian_standard_to_natural_params(mu_all, Sigma_all)
  #Q_joint <- rowSums(res_list$Lambda_all, dims = 2)
  Q_joint <- apply(res_list$Lambda_all, c(1,2), sum)
  r_joint <- apply(res_list$eta_all, 2, sum)
  vall <- f_log_norm_const_natural_gaussian_natural_params(r_joint, Q_joint)
  #browser()
  #i <- 1
  #res_inter <- f_gaussian_standard_to_natural_params(res_list$eta_all, res_list$Lambda_all)
  #normalise_product_norm(t(res_inter$eta_all), res_inter$Lambda_all)
  vsum <- sum(sapply(c(1:ssplits), function(i) f_log_norm_const_natural_gaussian_natural_params(res_list$eta_all[i,], res_list$Lambda_all[,,i])))
  return(vall-vsum)
}

if(F){
  v1 <- f_log_norm_const_natural_gaussian_standard_params(mu_all[1,], Sigma_all[,,1])
  v2 <- f_log_norm_const_natural_gaussian_standard_params(mu_all[2,], Sigma_all[,,2])
  v3 <- f_log_norm_const_natural_gaussian_standard_params(mu_all[3,], Sigma_all[,,3])

  #v1+v2+v3
  f_integral_product_gaussian(mu_all, Sigma_all)

}
