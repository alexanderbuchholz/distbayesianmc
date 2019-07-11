# iterate sampling

#library(MASS)
#library(BayesLogit)
#source("R/helperfunctions/f_helper.R")
#source("R/logisticreg/functions_logistic_reg.R")
#source("R/logisticreg/ep_logistic_regression.R")
#source("R/helperfunctions/f_dataset_loader_logistic_reg.R")


f_function_model_split_logistic <- function(ssplits, iter = 1, iseed = 42, dataset = "sim", typesplit = "random") {
  ### runs a logistic regression on the splitted data
  ###
  datalist <- f_dataset_loader(dataset)
  X <- datalist[["X"]]
  y <- datalist[["y"]]
  scale = 1
  datasplits <- f_pack_split_data(X, y, ssplits = ssplits, iseed = iseed, typesplit = typesplit)
  rm(X, y, datalist)
  datasplits <- f_prep_prior_logistic(datasplits, scale = scale)
  res <- f_logistic_splits(datasplits)
  return(res)
}



# load the data
f_function_model_split <- function(ssplits, iter = 1, iseed = 42, n_steps = 22000, burnin = 2000, dataset = "sim", typesplit = "random", returnres=F, scale=1) {
  datalist <- f_dataset_loader(dataset)
  X <- datalist[["X"]]
  y <- datalist[["y"]]
  datasplits <- f_pack_split_data(X, y, ssplits = ssplits, iseed = iseed, typesplit = typesplit)
  rm(X, y, datalist)
  #browser()
  #datasplits <- f_pack_all_data(X,y)
  datasplits <- f_prep_prior_logistic(datasplits, scale = scale)

  filename_individual_sim <- paste("individual_sim_", dataset, "_", ssplits, "_splits_rep_", iter, "_seed_", iseed, "_indrep", sep = "")
  filename_small <- paste("small_sim_", dataset, "_", ssplits, "_splits_rep_", iter, "_seed_", iseed, "_", typesplit, ".RData", sep = "")

  print("start sampling")
  set.seed(iter)
  res_chain <- f_parallel_sampler(datasplits, n_steps = n_steps, burnin = burnin, filename_individual_sim = filename_individual_sim, returnres = returnres)
  #browser()
  print("calc product normal")
  reslistisub <- logistic_isub_adapted(filename_individual_sim, ssplits)
  logsubposterior <- reslistisub$logIsub
  ESS <- reslistisub$ESS
  #alphasubpriorconstant <- ssplits*log(ssplits**(datasplits[[1]]$d/2)*( (det(2*pi*datasplits[[1]][["Bprior"]]))**0.5 )**(1-1/ssplits))
  print("calc sub prior constant")
  #browser()
  alphasubpriorconstant <- ssplits*f_alpha_sub(ssplits = ssplits, Vprior = datasplits[[1]]$Bprior/ssplits)

  logpostown <- 0
  vec_subposteriors <- rep(0, ssplits)
  vec_gaussianconsts <- rep(0, ssplits)
  vec_loglik <- rep(0, ssplits)
  mat_means <- matrix(0, ncol=datasplits[[1]]$d, nrow=ssplits)
  mat_cov <- array(0, dim=c(datasplits[[1]]$d, datasplits[[1]]$d, ssplits))
  for (i in 1:ssplits) {
    load(paste(filename_individual_sim, toString(i), ".RData", sep = ""))
    logpostown <- logpostown + srepchain$chib
    vec_subposteriors[i] <- srepchain$chib
    #browser()
    # approximation using the laplace metropolis estimator
    # calculate momentes
    cov_inter <- cov(t(srepchain$chainbeta))
    mean_post <- apply(srepchain$chainbeta, 1, mean)
    mat_means[i,] <- mean_post
    mat_cov[,,i] <- cov_inter

    # calculate log_ parts
    partdet <- 0.5*(determinant((cov_inter)))$modulus[1]
    part_lik <- f_loglik(mean_post, params = datasplits[[i]])
    vec_loglik[i] <- part_lik
    part_prior <- f_logprior(mean_post, params = datasplits[[i]])
    part_dim <- (datasplits[[i]]$d/2)*log(2*pi)
    approx_evidence <- part_lik+part_prior+partdet+part_dim
    vec_gaussianconsts[i] <- approx_evidence
    #browser()
    unlink(paste(filename_individual_sim, toString(i), ".RData", sep = ""))
  }
  #browser()
  # TODO:
  # 1. create aggregate mean, cov for natural params
  # 2. create full prior
  # put all the bits together
  part_all_normconst <- f_integral_product_gaussian(mat_means, mat_cov)
  #part_all_prior <-dmvnorm(x, mean= params[["bprior"]], sigma=params[["Bprior"]], log = TRUE)
  #Z_all <-


  print("start saving")
  res_small <- list()
  res_small[["alphasub"]] <- alphasubpriorconstant
  res_small[["logsubposterior"]] <- logsubposterior
  res_small[["logsubposteriorapprox"]] <- part_all_normconst
  res_small[["logpost"]] <- logpostown
  res_small[["normconstcombined"]] <- logpostown + logsubposterior + alphasubpriorconstant
  res_small[["normconstcombinedapprox"]] <- logpostown + part_all_normconst + alphasubpriorconstant
  res_small[["ESS"]] <- ESS
  res_small[["vec_logsubposteriors"]] <- vec_subposteriors
  res_small[["vec_condintegral"]] <- reslistisub$condintegral
  res_small[["vec_loggaussianconst"]] <- vec_gaussianconsts
  res_small[["mat_means"]] <- mat_means
  res_small[["mat_cov"]] <- mat_cov

  #res_chains[["normconstcombined"]] <- logpostown + logsubposterior + alphasubpriorconstant
  #res_chains[["logpost"]] <- logpostown
  #res_chains[["alphasub"]] <- alphasubpriorconstant
  #res_chains[["logsubposterior"]] <- logsubposterior


  #save(res_chains, file=filename_full)
  #browser()
  save(res_small, file = filename_small)
  if(returnres == F) return(res_small)
  else return(list(res_small, res_chain))
}

f_function_model_split_daniel <- function(ssplits, iter = 1, iseed = 42, n_steps = 22000, burnin = 2000, dataset = "sim", typesplit = "random") {
  datalist <- f_dataset_load_logistic_regression(dataset)
  X <- datalist[["X"]]
  y <- datalist[["y"]]
  scale = 1
  datasplits <- f_pack_split_data(X,y, ssplits = ssplits, iseed = iseed, typesplit = typesplit)
  rm(X, y, datalist)
  #datasplits <- f_pack_all_data(X,y)
  datasplits <- f_prep_prior_logistic(datasplits, scale = scale)
  set.seed(iter)

  srepchainlistdaniel <- list()
  for(i_s in c(1:ssplits)){
    ressample <- BayesLogit::logit(datasplits[[i_s]][["y"]], datasplits[[i_s]][["X"]], m0=datasplits[[i_s]][["bprior"]], P0=datasplits[[i_s]][["Bpriorinv"]], samp=n_steps, burn=burnin)

    #apply(ressample$beta, 2, mean)
    ressample <- augment_logitfit(ressample, datasplits[[i_s]][["bprior"]], datasplits[[i_s]][["Bprior"]])
    #logistic_isubres <- logistic_isub(list(ressample1, ressample2))
    betastar = apply(ressample$beta, 2, mean)

    chib <- chib_logitfit(ressample, datasplits[[i_s]][["bprior"]], datasplits[[i_s]][["Bprior"]], betastar)
    ressample[["chib"]] <- chib
    # lb <- rep(-Inf, datasplits[[i_s]]$d)
    # ub <- rep(Inf, datasplits[[i_s]]$d)
    # betaposterior <- ressample$beta#[burnin:n_steps, ]
    # colnames(betaposterior) <- colnames(datasplits[[i_s]][["X"]])
    # names(lb) <- names(ub) <- colnames(betaposterior)
    # bridge_result <- bridge_sampler(samples = as.matrix(betaposterior), data=NULL, log_posterior = f_logposterior, params =  datasplits[[i_s]], lb = lb, ub = ub, silent = TRUE)
    # ressample[["bridge"]] <- bridge_result$logml
    srepchainlistdaniel[[i_s]] <- ressample
  }

  reslistisub <- logistic_isub(srepchainlistdaniel)

  logsubposteriordaniel <- reslistisub$logIsub
  ESS <- reslistisub$ESS
  logpostdaniel <- 0
  vec_subposteriors <- rep(0, ssplits)
  i_counter = 1
  for(srep in srepchainlistdaniel) {
    logpostdaniel <- logpostdaniel + srep$chib
    vec_subposteriors[i_counter] <- srep$chib
    i_counter <- i_counter+1
  }
  alphasubpriorconstant <- ssplits * f_alpha_sub(ssplits = ssplits, Vprior = datasplits[[1]]$Bprior / ssplits)
  #logpostdaniel+logsubposteriordaniel+alphasubpriorconstant


  print("start saving")
  res_small <- list()
  res_small[["alphasub"]] <- alphasubpriorconstant
  res_small[["logsubposterior"]] <- logsubposteriordaniel
  res_small[["logpost"]] <- logpostdaniel
  res_small[["normconstcombined"]] <- logpostdaniel+logsubposteriordaniel+alphasubpriorconstant
  res_small[["ESS"]] <- ESS
  res_small[["vec_logsubposteriors"]] <- vec_subposteriors

  srepchainlistdaniel[["normconstcombined"]] <- logpostdaniel+logsubposteriordaniel+alphasubpriorconstant
  srepchainlistdaniel[["logpost"]] <- logpostdaniel
  srepchainlistdaniel[["alphasub"]] <- alphasubpriorconstant
  srepchainlistdaniel[["logsubposterior"]] <- logsubposteriordaniel

  filename_full <- paste("full_sim_", dataset, "_daniel_", ssplits, "_splits_rep_", iter, "_seed_", iseed, ".RData", sep = "")
  filename_small <- paste("small_sim_", dataset, "_daniel_", ssplits, "_splits_rep_", iter, "_seed_", iseed, "_", typesplit, ".RData", sep = "")
  #save(srepchainlistdaniel, file=filename_full)
  save(res_small, file=filename_small)
}
