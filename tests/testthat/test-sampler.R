context("test-sampler")

test_that("normalizing constant works", {

  data_pima = Pima.tr
  X = scale(data_pima[,1:7])
  ones = rep(1, nrow(data_pima))
  X = cbind(ones, X)
  y = (data_pima[,8] == "Yes")*1

  datasplits <- f_pack_split_data(X,y, ssplits = 1)
  datasplits <- f_prep_prior_logistic(datasplits)

  res_chain <- f_sample_chain(n_steps=10000, params=datasplits[[1]], burnin = 2000)

  negf_posterior <- function(x, params.=params) -f_loglik(x, params = params.)-f_logprior(x, params = params.)
  resoptim <- optim(datasplits[[1]]$bprior, negf_posterior, params=datasplits[[1]], control= c(maxit=10000), hessian=TRUE)

  # compare estimation of posterior mean/ map
  mse = mean((apply(res_chain$chainbeta, 1, mean)-resoptim$par)**2)
  expect_lt(mse, 0.1)

  covarproxy <- 2*diag(diag(solve(resoptim$hessian)))
  muproxy <- resoptim$par
  samples_is <- rmvnorm(n= 200000, mean= muproxy, sigma=covarproxy)
  logweights <- f_loglik(samples_is, params = datasplits[[1]])+f_logprior(samples_is, params = datasplits[[1]]) - dmvnorm(samples_is, mean=muproxy, sigma=covarproxy, log=TRUE)
  normconstis <- log(mean(exp(logweights)))

  ## Use Chib's method
  normconstchib <- f_logprior(muproxy, datasplits[[1]])-f_logconditional_mean(muproxy, res_chain, datasplits[[1]])+f_loglik(muproxy, datasplits[[1]])
  # compare normalizing constants
  expect_lt((normconstchib-normconstis)**2, 0.1)

  # using bridgesampling
  library(bridgesampling)
  lb <- rep(-Inf, datasplits[[1]][["d"]])
  ub <- rep(Inf, datasplits[[1]][["d"]])
  betaposterior <- data.frame(t(res_chain$chainbeta))
  colnames(betaposterior) <- colnames(datasplits[[1]][["X"]])
  names(lb) <- names(ub) <- colnames(betaposterior)
  bridge_result <- bridge_sampler(samples = as.matrix(betaposterior), data=NULL, log_posterior = f_logposterior, params =  datasplits[[1]], lb = lb, ub = ub, silent = TRUE)

  # compare normalizing constants
  expect_lt((normconstchib-bridge_result$logml)**2, 0.1)

})

test_that("stan+bridge sampling and data augmentation for the logistic reg work", {
  ssplits <- 5
  scale <- 1
  dataset <- "pima"
  iseed <- 42
  iter <- 1
  res_exact <- f_function_model_split(ssplits = ssplits, iter = iter, iseed = iseed, n_steps = 2000, burnin = 200, dataset = dataset, typesplit = "random", returnres = T, scale=scale)
  
  fileName <- "./stan_files/logistic_reg.stan"
  stan_code <- readChar(fileName, file.info(fileName)$size)
  mod <- stan_model(model_code = stan_code, auto_write = T)
  
  
  dataset_loaded <- f_dataset_loader(dataset)
  splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, P=0, ssplits=ssplits, iseed=42, typesplit="random")
  splitted_data <- f_prep_prior_logistic(splitted_data, scale = scale)
  
  res_approx <-  f_stan_sampling_splitted_data(mod, splitted_data, dataset = dataset, i_seed = iter, iter = iter, typesplit = "random", typeprior = "normal", nchain=2000)
  
  expect_lt((res_approx$normconstcombined-res_exact[[1]]$normconstcombined)**2, 2)
  expect_lt((res_approx$normconstcombined-res_exact[[1]]$normconstcombinedapprox)**2, 2)
})


test_that("subsampling works", {
  data_pima = Pima.tr
  X = scale(data_pima[,1:7])
  ones = rep(1, nrow(data_pima))
  X = cbind(ones, X)
  y = (data_pima[,8] == "Yes")*1
  ssplits= 5
  datasplits <- f_pack_split_data(X,y, ssplits = ssplits)
  datasplits <- f_prep_prior_logistic(datasplits)
  f_parallel_sampler(datasplits = datasplits)

  eplist <-f_ep_splits(datasplits)
  for(i in 1:ssplits){
    load(paste("tempsavesim", toString(i), ".RData", sep=""))
    expect_lt((eplist[[i]]$Z_exact-srepchain$chib)**2, 0.1)
    expect_lt((eplist[[i]]$Z_exact-eplist[[i]]$Z_ep)**2, 4) # this test is less precises
    #expect_lt((eplist[[i]]$Z_exact-srepchain$bridge)**2, 0.1)
    #expect_lt((srepchain$chib-srepchain$bridge)**2, 0.1)
    unlink(paste("tempsavesim", toString(i), ".RData", sep=""))
  }


})

test_that("stratification works", {
  ssplits <- 3
  res1 <- f_function_model_split(ssplits, iter = 1, iseed = 42, n_steps = 10000, burnin = 2000, dataset = "pima", typesplit = "random")
  res2 <- f_function_model_split(ssplits, iter = 1, iseed = 42, n_steps = 10000, burnin = 2000, dataset = "pima", typesplit = "strat_y")
  res3 <- f_function_model_split(ssplits, iter = 1, iseed = 42, n_steps = 10000, burnin = 2000, dataset = "pima", typesplit = "strat_y_cluster")
  expect_lt((res1[["normconstcombined"]] - res2[["normconstcombined"]])**2, 4) # this test is less precises
  expect_lt((res1[["normconstcombined"]] - res3[["normconstcombined"]])**2, 4) # this test is less precises
  unlink("small_sim_pima_3_splits_rep_1_seed_42_random.RData")
  unlink("small_sim_pima_3_splits_rep_1_seed_42_strat_y.RData")
  unlink("small_sim_pima_3_splits_rep_1_seed_42_strat_y_cluster.RData")
})

test_that("natural parameters work", {
  d <- 2
  ssplits <- 3
  mu_all <- rmvnorm(ssplits, mean = rep(0, d), sigma = diag(rep(1, d)))
  Sigma_all <- rWishart(n = ssplits, df = 5, Sigma = diag(rep(1, d)))

  # importance based sampling
  n <- 10**6
  is_samples <- rmvnorm(n, mean = mu_all[1,], sigma = Sigma_all[,,1])
  log_weights <- dmvnorm(is_samples, mean = mu_all[2,], sigma = Sigma_all[,,2], log = T)+dmvnorm(is_samples, mean = mu_all[3,], sigma = Sigma_all[,,3], log = T)
  is_approx1 = log_sum_exp(log_weights) - log(n)

  is_samples <- rmvnorm(n, mean = mu_all[2,], sigma = Sigma_all[,,2])
  log_weights <- dmvnorm(is_samples, mean = mu_all[1,], sigma = Sigma_all[,,1], log = T)+dmvnorm(is_samples, mean = mu_all[3,], sigma = Sigma_all[,,3], log = T)
  is_approx2 = log_sum_exp(log_weights) - log(n)

  is_samples <- rmvnorm(n, mean = mu_all[3,], sigma = Sigma_all[,,3])
  log_weights <- dmvnorm(is_samples, mean = mu_all[1,], sigma = Sigma_all[,,2], log = T)+dmvnorm(is_samples, mean = mu_all[3,], sigma = Sigma_all[,,3], log = T)
  is_approx3 = log_sum_exp(log_weights) - log(n)

  expect_lt((is_approx1-is_approx2)**2, 4) # this test is less precises
  expect_lt((is_approx3-is_approx2)**2, 4) # this test is less precises
  expect_lt((is_approx3-is_approx1)**2, 4) # this test is less precises

  # test natural parameter transformation
  res_list <- f_gaussian_standard_to_natural_params(mu_all, Sigma_all)
  res_list_inverse <- f_gaussian_standard_to_natural_params(res_list$eta_all, res_list$Lambda_all)
  expect_lt(sum((res_list_inverse$eta_all-mu_all)**2), 0.001)
  expect_lt(sum((res_list_inverse$Lambda_all-Sigma_all)**2), 0.001)

  # test log normalization constant
  const1 <- f_log_norm_const_natural_gaussian_standard_params(mu_all[1,], Sigma_all[,,1])
  const2 <- f_log_norm_const_natural_gaussian_natural_params(res_list$eta_all[1,], res_list$Lambda_all[,,1])
  expect_lt((const1-const2)*2, 0.00001)

  # test density functions
  n = 10**1
  normal_is <- rmvnorm(n, mean=mu_all[1,], sigma = Sigma_all[,,1])

  diff_density <- f_exponential_form_normal(normal_is[1,], mu_all[1,], Sigma_all[,,1])-
    f_log_norm_const_natural_gaussian_standard_params(mu_all[1,], Sigma_all[,,1])-
    dmvnorm(normal_is[1,], mean=mu_all[1,], sigma = Sigma_all[,,1], log = T)
  expect_lt((diff_density)**2, 0.000001)

  # test product of integral
  true_const <- f_integral_product_gaussian(mu_all, Sigma_all)
  approx_const <- (is_approx1+is_approx2+is_approx3)/3
  expect_lt((true_const-approx_const)**2, 1) # as we use IS we have to allow for some MC error


})
