# params simulation linear model with laplace prior

# specifics for stan
fileName <- "~/R_programming/distbayesianmc/stan_files/sparse_linear_reg_laplace.stan"

# single simulation

scale <- 1

#dataset <- "higgs1"
#iseed <- 42
#iter <- 1
#ssplits <- 5

typesplit <-  "random"
nobs <- 4*10**3
nchain <- 2000
typeprior <- "laplace_normal"
highcorr <- TRUE

Pparams <- 2

iters <- 10
#vec_datasets <- c("sparse_reg_1", "sparse_reg_2")
vec_datasets <- c("hla_ultra_small1", "hla_ultra_small2")
vec_types_splits <- c(typesplit)

# multiple simulations
vec_splits <- c(1,5,10,20,50)#1,30,50)#,10,20,50)
