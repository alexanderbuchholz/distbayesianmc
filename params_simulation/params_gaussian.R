# params simulation linear model with laplace prior

# specifics for stan
fileName <- "~/R_programming/distbayesianmc/stan_files/linear_reg_gaussian.stan"

# single simulation

scale <- 1

#dataset <- "higgs1"
#iseed <- 42
#iter <- 1
#ssplits <- 5

typesplit <-  "random"
nobs <- 10**4
nchain <- 10000
typeprior <- "normal"
highcorr <- TRUE

Pparams <- 2

iters <- 20
#vec_datasets <- c("sparse_reg_1", "sparse_reg_2")
#vec_datasets <- c("hla1", "hla2")
vec_datasets <- c("gaussian1", "gaussian2", "gaussian3","gaussian4", "gaussian5", "gaussian6") # 
vec_types_splits <- c(typesplit)

# multiple simulations
vec_splits <- c(1,2,5,10,20,50)#1,30,50)#,10,20,50)
