# transform the datasplit to a stan model

f_stan_sampling_splitted_data_logistic <- function(mod, splitted_data, dataset, i_seed, iter, typesplit, epapprox=T, bridgepack=F, typeprior="normal", nchain=2000, file_identifier=""){
  ssplits <- length(splitted_data)
  res_stan_sampling <- list()
  
  
  mat_means <- matrix(0, ncol=splitted_data[[1]]$d, nrow=ssplits)
  mat_cov <- array(0, dim=c(splitted_data[[1]]$d, splitted_data[[1]]$d, ssplits))
  
  alphasubpriorconstant <- ssplits*f_alpha_sub(ssplits = ssplits, Vprior = splitted_data[[1]]$Bprior/ssplits, typeprior = typeprior)
  vec_logsubpost <- rep(0, ssplits)
  vec_logsubpost_ep <- rep(0, ssplits)
  beta_list <- list()
  for(i_s in 1:ssplits){
    res_stan_sampling[[i_s]] <- list()
    dat <- list(N        = splitted_data[[i_s]]$n,
                P        = splitted_data[[i_s]]$d,
                y    = splitted_data[[i_s]]$y,
                x = splitted_data[[i_s]]$X,
                sigma = splitted_data[[i_s]]$scale
    )
    #resStan <- stan(model_code = stan_code, data = dat,
    #                chains = 1, iter = 3000, warmup = 500, thin = 1, auto_write=T)
    
    resStan <- rstan::sampling(mod, data = dat, chains = 1, iter = nchain, warmup = nchain*0.2, thin = 1, seed = i_seed)
    interres <- rstan::extract(resStan, pars="beta", permuted=F)
    
    if(epapprox){
      epres <- EPlogit(Y=splitted_data[[i_s]]$y, X = splitted_data[[i_s]]$X, s = splitted_data[[i_s]]$scale)
      vec_logsubpost_ep[i_s] <- epres$Z 
    }
    
    if(bridgepack){
      bridge_results <- bridge_sampler(samples= resStan, silent = F, stanfit_model = resStan, method = "warp3") 
    }
    else{
      lb <- rep(-Inf, splitted_data[[i_s]][["d"]])
      ub <- rep(Inf, splitted_data[[i_s]][["d"]])
      betaposterior <- array(interres, c(dim(interres)[1], dim(interres)[3]))
      colnames(betaposterior) <- colnames(splitted_data[[i_s]][["X"]])
      names(lb) <- names(ub) <- colnames(betaposterior)
      #browser()
      bridge_results <- bridge_sampler(samples = as.matrix(betaposterior), data=NULL, log_posterior = f_logposterior, params =  splitted_data[[i_s]], lb = lb, ub = ub, silent = TRUE)
    }
    
    
    res_stan_sampling[[i_s]][["subposteriorevidence"]] <- bridge_results$logml
    res_stan_sampling[[i_s]][["subposteriormean"]] <- apply(interres, 3, mean)
    vec_logsubpost[i_s] <- bridge_results$logml
    res_stan_sampling[[i_s]][["subposteriorvariance"]] <- cov(array(interres, c(dim(interres)[1], dim(interres)[3])))
    mat_means[i_s,] <- res_stan_sampling[[i_s]][["subposteriormean"]]
    mat_cov[,,i_s] <- res_stan_sampling[[i_s]][["subposteriorvariance"]]
    betasamples <- array(interres, c(dim(interres)[1], dim(interres)[3]))
    beta_list[[i_s]] <- betasamples
  }
  part_all_normconst <- f_integral_product_gaussian(mat_means, mat_cov)
  res_small <- list()
  logsubpost <- sum(vec_logsubpost)
  res_small[["alphasub"]] <- alphasubpriorconstant
  res_small[["logsubposteriorapprox"]] <- part_all_normconst
  res_small[["logpost"]] <- logsubpost
  res_small[["vec_logsubpost"]] <- vec_logsubpost
  res_small[["vec_logsubpost_ep"]] <- vec_logsubpost_ep
  res_small[["normconstcombined"]] <- logsubpost + part_all_normconst + alphasubpriorconstant
  res_small[["normconstcombined_type"]] <- "approx"
  res_small[["mat_means"]] <- mat_means
  res_small[["mat_cov"]] <- mat_cov
  res_small[["betasamples"]] <- beta_list
  
  filename_small <- paste("small_sim_stan_", dataset, "_", ssplits, "_splits_rep_", iter, "_seed_", i_seed, "_", typesplit, file_identifier, ".RData", sep = "") 
  save(res_small, file = filename_small)
  return(res_small)
}



f_stan_sampling_splitted_data <- function(mod, splitted_data, dataset, i_seed, iter, typesplit, bridgepack=T, typeprior="laplace_normal", nchain=2000, file_identifier=""){
  ssplits <- length(splitted_data)
  res_stan_sampling <- list()
  
  
  mat_means <- matrix(0, ncol=splitted_data[[1]]$Pparams, nrow=ssplits)
  mat_cov <- array(0, dim=c(splitted_data[[1]]$Pparams, splitted_data[[1]]$Pparams, ssplits))
  
  #browser()
  alphasubpriorconstant <- ssplits*f_alpha_sub(ssplits = ssplits, Vprior = splitted_data[[1]]$Bprior/ssplits, typeprior = typeprior)
  vec_logsubpost <- rep(0, ssplits)
  vec_logsubpost_is <- rep(0, ssplits)
  vec_logsubpost_ep <- rep(0, ssplits)
  beta_list <- list()
  for(i_s in 1:ssplits){
    res_stan_sampling[[i_s]] <- list()
    dat <- list(N        = splitted_data[[i_s]]$n,
                P        = splitted_data[[i_s]]$d,
                y    = splitted_data[[i_s]]$y,
                x = splitted_data[[i_s]]$X,
                sigma = splitted_data[[i_s]]$scale
    )
    #resStan <- stan(model_code = stan_code, data = dat,
    #                chains = 1, iter = 3000, warmup = 500, thin = 1, auto_write=T)
    res_stan_sampling
    resStan <- rstan::sampling(mod, data = dat, chains = 1, iter = nchain, warmup = nchain*0.2, thin = 1, seed = i_seed)
    interres <- rstan::extract(resStan, pars="beta", permuted=F)
    
    
    if(bridgepack){
      bridge_results <- bridge_sampler(samples= resStan, silent = F, stanfit_model = resStan, method = "warp3") 
    }
    else{
      lb <- rep(-Inf, splitted_data[[i_s]][["d"]])
      ub <- rep(Inf, splitted_data[[i_s]][["d"]])
      betaposterior <- array(interres, c(dim(interres)[1], dim(interres)[3]))
      colnames(betaposterior) <- colnames(splitted_data[[i_s]][["X"]])
      names(lb) <- names(ub) <- colnames(betaposterior)
      #browser()
      bridge_results <- bridge_sampler(samples = as.matrix(betaposterior), data=NULL, log_posterior = f_logposterior, params =  splitted_data[[i_s]], lb = lb, ub = ub, silent = TRUE)
    }
    
    
    res_stan_sampling[[i_s]][["subposteriorevidence"]] <- bridge_results$logml
    res_stan_sampling[[i_s]][["subposteriormean"]] <- apply(interres, 3, mean)
    vec_logsubpost[i_s] <- bridge_results$logml
    res_stan_sampling[[i_s]][["subposteriorvariance"]] <- cov(array(interres, c(dim(interres)[1], dim(interres)[3])))
    mat_means[i_s,] <- res_stan_sampling[[i_s]][["subposteriormean"]]
    mat_cov[,,i_s] <- res_stan_sampling[[i_s]][["subposteriorvariance"]]
    betasamples <- array(interres, c(dim(interres)[1], dim(interres)[3]))
    beta_list[[i_s]] <- betasamples
    #resclust <- Mclust(betasamples, G=2)
    #restest <- sapply(1:splitted_data[[i_s]]$Pparams,  function(i) shapiro.test(betasamples[,i])$p.value)
    #browser()
    ### IS sampling for normalizing constant, for testing that our code does the right thing
    if(F){
      beta_proposal <- mvrnorm(n = 20, mu = mat_means[i_s,], Sigma = mat_cov[,,i_s])
      lw_prop <- dmvnorm(x=beta_proposal, mean = mat_means[i_s,], sigma = mat_cov[,,i_s], log = T)
      lwprior <- f_logprior(beta_proposal, splitted_data[[i_s]])
      
      n_prop <- dim(beta_proposal)[1]
      lwlik <- rep(0, n_prop)
      for(i_n in 1:n_prop){
        lwlik[i_n] <- sum(dmvnorm(x=splitted_data[[i_s]]$y, mean = ((beta_proposal[i_n,1:10, drop=F]) %*% t(splitted_data[[i_s]]$X))+beta_proposal[i_n,11], log = TRUE))
      }
      vec_logsubpost_is[i_s] <- log_sum_exp(lwprior+lwlik-lw_prop)-log(n_prop)
      #browser()
    }
  }
  part_all_normconst <- f_integral_product_gaussian(mat_means, mat_cov)
  res_small <- list()
  logsubpost <- sum(vec_logsubpost)
  res_small[["alphasub"]] <- alphasubpriorconstant
  res_small[["logsubposteriorapprox"]] <- part_all_normconst
  res_small[["logpost"]] <- logsubpost
  res_small[["vec_logsubpost"]] <- vec_logsubpost
  res_small[["vec_logsubpost_ep"]] <- vec_logsubpost_ep
  res_small[["vec_logsubpost_is"]] <- vec_logsubpost_is
  res_small[["normconstcombined"]] <- logsubpost + part_all_normconst + alphasubpriorconstant
  res_small[["normconstcombined_type"]] <- "approx"
  res_small[["mat_means"]] <- mat_means
  res_small[["mat_cov"]] <- mat_cov
  res_small[["betasamples"]] <- beta_list
  
  filename_small <- paste("small_sim_stan_", dataset, "_", ssplits, "_splits_rep_", iter, "_seed_", i_seed, "_", typesplit, file_identifier, ".RData", sep = "") 
  save(res_small, file = filename_small)
  return(res_small)
}


f_stan_sampling_splitted_data_median <- function(mod, splitted_data, dataset, i_seed, iter, typesplit, bridgepack=T, typeprior="laplace_normal", nchain=2000, file_identifier=""){
  ssplits <- length(splitted_data)
  res_stan_sampling <- list()
  
  
  mat_means <- matrix(0, ncol=splitted_data[[1]]$Pparams, nrow=ssplits)
  mat_cov <- array(0, dim=c(splitted_data[[1]]$Pparams, splitted_data[[1]]$Pparams, ssplits))
  
  #browser()
  alphasubpriorconstant <- ssplits*f_alpha_sub(ssplits = ssplits, Vprior = splitted_data[[1]]$Bprior/ssplits, typeprior = typeprior)
  vec_logsubpost <- rep(0, ssplits)
  vec_logsubpost_is <- rep(0, ssplits)
  vec_logsubpost_ep <- rep(0, ssplits)
  beta_list <- list()
  for(i_s in 1:ssplits){
    res_stan_sampling[[i_s]] <- list()
    dat <- list(N        = splitted_data[[i_s]]$n,
                P        = splitted_data[[i_s]]$d,
                y    = splitted_data[[i_s]]$y,
                x = splitted_data[[i_s]]$X,
                sigma = splitted_data[[i_s]]$scale,
                R = splitted_data[[i_s]]$R
    )
    #resStan <- stan(model_code = stan_code, data = dat,
    #                chains = 1, iter = 3000, warmup = 500, thin = 1, auto_write=T)
    res_stan_sampling
    resStan <- rstan::sampling(mod, data = dat, chains = 1, iter = nchain, warmup = nchain*0.2, thin = 1, seed = i_seed)
    interres <- rstan::extract(resStan, pars="beta", permuted=F)
    
    
    if(bridgepack){
      bridge_results <- bridge_sampler(samples= resStan, silent = F, stanfit_model = resStan, method = "warp3") 
    }
    else{
      lb <- rep(-Inf, splitted_data[[i_s]][["d"]])
      ub <- rep(Inf, splitted_data[[i_s]][["d"]])
      betaposterior <- array(interres, c(dim(interres)[1], dim(interres)[3]))
      colnames(betaposterior) <- colnames(splitted_data[[i_s]][["X"]])
      names(lb) <- names(ub) <- colnames(betaposterior)
      #browser()
      bridge_results <- bridge_sampler(samples = as.matrix(betaposterior), data=NULL, log_posterior = f_logposterior, params =  splitted_data[[i_s]], lb = lb, ub = ub, silent = TRUE)
    }
    
    
    res_stan_sampling[[i_s]][["subposteriorevidence"]] <- bridge_results$logml
    res_stan_sampling[[i_s]][["subposteriormean"]] <- apply(interres, 3, mean)
    vec_logsubpost[i_s] <- bridge_results$logml
    res_stan_sampling[[i_s]][["subposteriorvariance"]] <- cov(array(interres, c(dim(interres)[1], dim(interres)[3])))
    mat_means[i_s,] <- res_stan_sampling[[i_s]][["subposteriormean"]]
    mat_cov[,,i_s] <- res_stan_sampling[[i_s]][["subposteriorvariance"]]
    betasamples <- array(interres, c(dim(interres)[1], dim(interres)[3]))
    beta_list[[i_s]] <- betasamples
    #resclust <- Mclust(betasamples, G=2)
    #restest <- sapply(1:splitted_data[[i_s]]$Pparams,  function(i) shapiro.test(betasamples[,i])$p.value)
    #browser()
    ### IS sampling for normalizing constant, for testing that our code does the right thing
  }
  part_all_normconst <- f_integral_product_gaussian(mat_means, mat_cov)
  res_small <- list()
  logsubpost <- sum(vec_logsubpost)
  res_small[["alphasub"]] <- alphasubpriorconstant
  res_small[["logsubposteriorapprox"]] <- part_all_normconst
  res_small[["logpost"]] <- logsubpost
  res_small[["vec_logsubpost"]] <- vec_logsubpost
  res_small[["vec_logsubpost_ep"]] <- vec_logsubpost_ep
  res_small[["vec_logsubpost_is"]] <- vec_logsubpost_is
  res_small[["normconstcombined"]] <- logsubpost + part_all_normconst + alphasubpriorconstant
  res_small[["normconstcombined_type"]] <- "approx"
  res_small[["mat_means"]] <- mat_means
  res_small[["mat_cov"]] <- mat_cov
  res_small[["betasamples"]] <- beta_list
  
  filename_small <- paste("small_sim_stan_median_", dataset, "_", ssplits, "_splits_rep_", iter, "_seed_", i_seed, "_", typesplit, file_identifier, ".RData", sep = "") 
  save(res_small, file = filename_small)
  return(res_small)
}


f_stan_sampling_single_split <- function(mod, single_split, dataset, i_seed, iter, typesplit, bridgepack=T, typeprior="laplace_normal", nchain=2000, file_identifier=""){
  # run the stan sampler on a single split, needed for the higgs dataset
  start_time = Sys.time()
  
  res_stan_sampling <- list()
  
  ssplits <- single_split$ssplits
  #browser()
  alphasubpriorconstant <- ssplits*f_alpha_sub(ssplits = ssplits, Vprior = single_split$Bprior/ssplits, typeprior = typeprior)
  
  dat <- list(N        = single_split$n,
              P        = single_split$d,
              y    = single_split$y,
              x = single_split$X,
              sigma = single_split$scale
  )
  
  resStan <- rstan::sampling(mod, data = dat, chains = 1, iter = nchain, warmup = nchain*0.2, thin = 1, seed = i_seed)
  interres <- rstan::extract(resStan, pars="beta", permuted=F)
  print(tempdir())
  
  
  bridge_results <- bridge_sampler(samples= resStan, silent = F, stanfit_model = resStan, method = "warp3") 
  
  res_stan_sampling[["subposteriorevidence"]] <- bridge_results$logml
  res_stan_sampling[["subposteriormean"]] <- apply(interres, 3, mean)
  logsubpost <- bridge_results$logml
  res_stan_sampling[["subposteriorvariance"]] <- cov(array(interres, c(dim(interres)[1], dim(interres)[3])))
  betasamples <- array(interres, c(dim(interres)[1], dim(interres)[3]))

  res_small <- list()
  
  res_small[["alphasub"]] <- alphasubpriorconstant
  res_small[["logsubpost"]] <- logsubpost
  res_small[["normconstcombined_type"]] <- "approx"
  res_small[["mat_means"]] <- res_stan_sampling[["subposteriormean"]]
  res_small[["mat_cov"]] <- res_stan_sampling[["subposteriorvariance"]]
  res_small[["betasamples"]] <- betasamples
  
  end_time = Sys.time()
  timediff <- end_time-start_time
  res_small[["runtime"]] <- as.numeric(timediff)
  
  filename_small <- paste("small_sim_stan_ind_", dataset, "_", ssplits, "_splits_i_split_", single_split$isplit, "_rep_", iter, "_seed_", i_seed, "_", typesplit, file_identifier, ".RData", sep = "") 
  save(res_small, file = filename_small)
  return(res_small)
}
