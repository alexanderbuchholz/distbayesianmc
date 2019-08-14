# params simuatlion logistic

# specifics for stan
fileName <- "~/R_programming/distbayesianmc/stan_files/logistic_reg.stan"

# single simulation

scale <- 1

#dataset <- "higgs1"
#iseed <- 42
#iter <- 1
#ssplits <- 5

typesplit <-  "random"
nobs <- 10**3
nchain <- 3000
typeprior <- "normal"
Pparams <- 0
iters <- 20

#vec_datasets <- c("higgs1_full", "higgs2_full")
vec_datasets <- c("higgs1_large", "higgs2_large")#, "higgs1_full", "higgs2_full")
#vec_datasets <- c("pima")
#vec_types_splits <- c(typesplit, "strat_y_cluster")
vec_types_splits <- c(typesplit)

# multiple simulations
#vec_splits <- c(1,2,3,5,10,20)
#vec_splits <- c(100,500,1000)#,100)
#vec_splits <- c(2,3,5,10)#,20)
#vec_splits <- c(1)
