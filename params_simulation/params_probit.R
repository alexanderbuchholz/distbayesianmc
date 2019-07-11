# params simuatlion logistic

# specifics for stan
fileName <- "~/R_programming/distbayesianmc/stan_files/probit_reg.stan"

# single simulation

scale <- 1

dataset <- "pima"
iseed <- 42
iter <- 1
ssplits <- 5

typesplit <-  "random"
nobs <- 10**3
nchain <- 2000
typeprior <- "normal"

iters <- 10
vec_datasets <- c(dataset)
#vec_datasets <- c(dataset)
vec_types_splits <- c(typesplit)

# multiple simulations
vec_splits <- c(1,2,3,5,10)#,20)
