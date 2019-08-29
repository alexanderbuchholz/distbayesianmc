# help files for splitting

# helper for simulations on the cluster
f_i_iter <- function(sim_id, ssplits){
  # returns the iteration number (in thousands or hundreds)
  if(sim_id %% ssplits == 0){
    return((sim_id %/% ssplits))
  }
  else{
    return( (sim_id %/% ssplits) + 1)
  }
}
f_split_number <- function(sim_id, ssplits) {
  # function that returns the ssplit number
  if (sim_id <= ssplits){
    return(sim_id)
  }
  else {
    if((sim_id %% ssplits)==0){
      return(ssplits)
    }
    else{
      return((sim_id %% ssplits))
    }
    
  }
}


f_pack_split_data <- function(X, y, d=NULL, Pparams=0, ssplits=4, iseed=42, typesplit="random", dataset = ""){
  # returns a list with the splitted data
  set.seed(iseed)
  splitslist <- list()
  if(typesplit == "random") {
    print("random sampling")
    selector <- sample(1:ssplits, nrow(X), replace=T)
    #browser()
    for(i_s in c(1:ssplits)) {
      splitslist[[i_s]] <- list()
      splitslist[[i_s]][["X"]] <- X[selector == i_s,]
      splitslist[[i_s]][["y"]] <- y[selector == i_s]
      splitslist[[i_s]][["kappa"]] <- y[selector == i_s]-0.5
      if(!is.null(d)){
        splitslist[[i_s]][["d"]] <- d
      }
      else{
        splitslist[[i_s]][["d"]] <- dim(X[selector == i_s,,drop = F])[2]
      }
      splitslist[[i_s]][["Pparams"]] <- splitslist[[i_s]][["d"]]+Pparams
      splitslist[[i_s]][["n"]] <- dim(X[selector == i_s,,drop = F])[1]
      splitslist[[i_s]][["dataset"]] <- dataset
    }
  }
  else if(typesplit == "strat_y") {
    print("stratified sampling on y")
    print("try corrected version")
    df <- as.data.frame(cbind(X, y))
    numvars <- dim(X)[2]
    #browser()
    for(i_s in c(1:ssplits)) {
      if(ssplits == 1){
        out <- df
      }

      else{
        frac_split = (1/(ssplits - i_s + 1))
        if(frac_split == 1){
          out <- df
          assert("frac split = 1 although not on last split", {(ssplits==i_s)})
          #browser()
        }
        else{
          outlist <- splitstackshape::stratified(df, group="y", size = frac_split, bothSets = T)
          out <- outlist[["SAMP1"]]
          rm(df)
          df <- outlist[["SAMP2"]]
          #out <- df %>%
          #  group_by(y) %>%
          #  sample_frac(frac_split) %>% ungroup()

          #out <- out[sample(nrow(out)),]
        }

      }

      #out <- df %>%
      #  group_by(y) %>%
      #  sample_frac(1/(ssplits - i_s + 1)) %>% ungroup()
      #browser()
      mat <- data.matrix(out)
      splitslist[[i_s]] <- list()
      splitslist[[i_s]][["X"]] <- mat[,1:numvars]
      splitslist[[i_s]][["y"]] <- mat[,numvars + 1]
      splitslist[[i_s]][["kappa"]] <- mat[,numvars + 1] - 0.5
      splitslist[[i_s]][["d"]] <- numvars
      splitslist[[i_s]][["n"]] <- length(mat[,numvars + 1])
      splitslist[[i_s]][["Pparams"]] <- splitslist[[i_s]][["d"]]+Pparams
      splitslist[[i_s]][["dataset"]] <- dataset
      #df <- setdiff(df, out)
    }
  }
  else if(typesplit == "strat_y_cluster") {
    print("stratified sampling on y and clusters")
    resclustering <- kmeans(X, 20, iter.max = 50, nstart = 10)
    vec_clusters <- resclustering$cluster
    df <- as.data.frame(cbind(X, y, vec_clusters))
    #browser()
    numvars <- dim(X)[2]
    for(i_s in c(1:ssplits)) {
      if(ssplits == 1){
        out <- df
      }
      else{
        frac_split = 1/(ssplits - i_s + 1)
        if(frac_split == 1){
          out <- df
          assert("frac split = 1 although not on last split", {(ssplits==i_s)})
        }
        else{
          outlist <- splitstackshape::stratified(df, group=c("y", "vec_clusters"), size = frac_split, bothSets = T)
          out <- outlist[["SAMP1"]]
          rm(df)
          df <- outlist[["SAMP2"]]

          #out <- df %>%
          #  group_by(y, vec_clusters) %>%
          #  sample_frac(frac_split) %>% ungroup()

        }

      }



      mat <- data.matrix(out)
      splitslist[[i_s]] <- list()
      splitslist[[i_s]][["X"]] <- mat[,1:numvars]
      splitslist[[i_s]][["y"]] <- mat[,numvars + 1]
      splitslist[[i_s]][["kappa"]] <- mat[,numvars + 1]-0.5
      splitslist[[i_s]][["d"]] <- numvars
      splitslist[[i_s]][["n"]] <- length(mat[,numvars + 1])
      splitslist[[i_s]][["Pparams"]] <- splitslist[[i_s]][["d"]]+Pparams
      splitslist[[i_s]][["dataset"]] <- dataset
      #df <- setdiff(df, out)
    }
  }
  else if(typesplit == "strat_y_cluster_qmc") {
    print("stratified sampling on y and qmc based clusters")
    nobs <-  dim(X)[1]
    numvars <- dim(X)[2]
    splitsize  <-  floor(0.1*nobs/ssplits)
    covmat <- cov(X[,2:numvars])
    qmc_normal <- qnorm(randtoolbox::runif.halton(floor(nobs/ssplits), dim = numvars-1))
    #qmc_normal <- qnorm(randtoolbox::runif.sobol(floor(nobs/ssplits), dim = dim-1, scrambling = 1))
    qmc_normal <- qmc_normal%*%chol(covmat)
    distmat = pairwiseDistC(X[,2:numvars],qmc_normal)
    cluster_assignments = apply(distmat, 2, which.min)
    #browser()
    df <- as.data.frame(cbind(X, y, cluster_assignments))
    #

    for(i_s in c(1:ssplits)) {
      if(ssplits == 1){
        out <- df
      }
      else{
        frac_split = 1/(ssplits - i_s + 1)
        if(frac_split == 1){
          out <- df
          assert("frac split = 1 although not on last split", {(ssplits==i_s)})
        }
        else{
          outlist <- splitstackshape::stratified(df, group=c("y", "cluster_assignments"), size = frac_split, bothSets = T)
          out <- outlist[["SAMP1"]]
          rm(df)
          df <- outlist[["SAMP2"]]

          #out <- df %>%
          #  group_by(y, vec_clusters) %>%
          #  sample_frac(frac_split) %>% ungroup()

        }

      }



      mat <- data.matrix(out)
      splitslist[[i_s]] <- list()
      splitslist[[i_s]][["X"]] <- mat[,1:numvars]
      splitslist[[i_s]][["y"]] <- mat[,numvars + 1]
      splitslist[[i_s]][["kappa"]] <- mat[,numvars + 1]-0.5
      splitslist[[i_s]][["d"]] <- numvars
      splitslist[[i_s]][["n"]] <- length(mat[,numvars + 1])
      splitslist[[i_s]][["dataset"]] <- dataset
      #df <- setdiff(df, out)
    }
  }
  else if(typesplit == "strat_y_x") {
    print("stratified sampling on y and x")


    numvars <- dim(X)[2]
    median_vec <- apply(X, 2, median)
    compareX <- X[,2:numvars] < median_vec[2:numvars]
    namescompX <- paste("strat_", colnames(compareX), sep="")
    colnames(compareX) <- namescompX
    df <- as.data.frame(cbind(X, y, compareX))
    for(i_s in c(1:ssplits)) {
      #browser()
      if(ssplits == 1){
        out <- df
      }
      else{
        frac_split = 1/(ssplits - i_s + 1)
        if(frac_split == 1) {
          out <- df
          assert("frac split = 1 although not on last split", {(ssplits==i_s)})
        }
        else{
          outlist <- splitstackshape::stratified(df, group=c("y", namescompX), size = frac_split, bothSets = T)
          out <- outlist[["SAMP1"]]
          rm(df)
          df <- outlist[["SAMP2"]]



          #out <- df %>%
          #  group_by(.dots=c("y",namescompX)) %>%
          #  sample_frac(frac_split) %>% ungroup()

        }

      }


      mat <- data.matrix(out)
      splitslist[[i_s]] <- list()
      splitslist[[i_s]][["X"]] <- mat[,1:numvars]
      splitslist[[i_s]][["y"]] <- mat[,numvars + 1]
      splitslist[[i_s]][["kappa"]] <- mat[,numvars + 1]-0.5
      splitslist[[i_s]][["d"]] <- numvars
      splitslist[[i_s]][["n"]] <- length(mat[,numvars + 1])
      splitslist[[i_s]][["dataset"]] <- dataset
      #df <- setdiff(df, out)
    }
  }
  else if(typesplit == "disk_sampling") {
    print("disk sampling")
    nobs <-  dim(X)[1]
    d <-  dim(X)[2]

    covmat <- cov(X[,2:d])
    X_unif <- pnorm(X[,2:d] %*%  (solve(covmat)))
    r0 <- (0.75*ssplits/nobs*gamma(d/2 + 1)/(pi**(d/2)))**(1/d)
    ind_split <- cbind(rep(0, nobs), 1:nobs)
    ind_split[1:ssplits,1] <- c(1:ssplits)
    while(sum(ind_split[,1] == 0) > 0){
      # loop over the other indices, select an observation that is still available
      selector <- sample(ind_split[ind_split[,1] == 0, 2], 1)
      flag_removed <- FALSE
      min_inter_distance <- rep(0, ssplits)
      for (jrem in 1:ssplits){
        # check if fits into bucket
        distances <- sqrt(pairwiseDistC(t(as.matrix(X_unif[selector,])), matrix(X_unif[ind_split[,1] == jrem,], ncol = d-1)))
        min_inter_distance[jrem] <- min(distances)
        if ( all(distances > r0)){
          if(sum(ind_split[,1] == jrem) < ceiling(nobs/ssplits)){
            ind_split[selector,1] <- jrem
            flag_removed <- TRUE
          }

        }
        if((jrem == ssplits) & !(flag_removed)){
          ind_split[selector,1] <- which.max(min_inter_distance)
        }

      }
    }

    for(i_s in c(1:ssplits)) {
      selector <- ind_split[,1]
      splitslist[[i_s]] <- list()
      splitslist[[i_s]][["X"]] <- X[selector == i_s,]
      splitslist[[i_s]][["y"]] <- y[selector == i_s]
      splitslist[[i_s]][["kappa"]] <- y[selector == i_s]-0.5
      splitslist[[i_s]][["d"]] <- dim(X[selector == i_s,])[2]
      splitslist[[i_s]][["n"]] <- dim(X[selector == i_s,])[1]
      splitslist[[i_s]][["dataset"]] <- dataset
    }
  }
  return(splitslist)
}
