# params simuatlion logistic

# specifics for stan
fileName <- "./stan_files/logistic_reg.stan"

# single simulation
ssplits <- 5
scale <- 1
dataset <- "higgs1"
iseed <- 42
iter <- 1
typesplit <-  "random"
nobs <- 10**3
nchain <- 2000
typeprior <- "normal"
iters <- 10
vec_datasets <- c("higgs1", "higgs2")
#vec_datasets <- c(dataset)
vec_types_splits <- c(typesplit)

# multiple simulations
vec_splits <- c(1,2,3,5,10,20)
