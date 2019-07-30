# R2BGLiMS
library(R2BGLiMS)
library(distbayesianmc)

dataset_loaded <- f_dataset_loader("pima")
pvars <- dim(dataset_loaded$X)[2]
#pvars <- 4
#splitted_data <- f_pack_split_data(dataset_loaded$X, dataset_loaded$y, ssplits=ssplits, iseed=iter, typesplit="random")
X <-  dataset_loaded$X[,2:pvars]
predictors <- colnames(X)

data <- cbind(X, dataset_loaded$y)
colnames(data)[pvars] <- "outcome"

pima.results <- R2BGLiMS(
      likelihood="Logistic",
      data=data,
      outcome.var="outcome",
      g.prior = FALSE,
      model.space.priors=list(
          "a"=1, "b"=1, "Variables"=predictors),
      n.mil=5,
      seed=2,
      thinning.interval = 1,
      beta.prior.partitions = list(list("Variables" = predictors, "UniformA"=0.9999, "UniformB"=1)), # fix the hyper prior
      extra.arguments = list("Adaption_Iterations" = 2e5, 
                             "AlphaPriorSd" = 1,
                             "GaussianResidualPriorFamily" = 1,
                             "GaussianResidualPrior_UnifArg1" = 0.9999,
                             "GaussianResidualPrior_UnifArg2" = 1.
                             )
    )
topmodel_res <- TopModels(pima.results)
ManhattanPlot(pima.results)
pima.results@mcmc.output

#pima.results@posterior.summary.table

# count_full_collumns <- function(x){
#   all(!(x  == 0))
# }
# 
# count_collumns_not_bp_skin <- function(x){
#   selector_small <-  !(colnames(x) %in% c("bp", "skin", "age"))
#   return(all(!(x[selector_small]  == 0)) & all(x[!selector_small]  == 0))
# }


#count_full_collumns(pima.results@mcmc.output[9,])
#count_collumns_not_bp_skin(pima.results@mcmc.output[235,])

#mean(apply(pima.results@mcmc.output, 1, count_full_collumns))
#mean(apply(pima.results@mcmc.output, 1, count_collumns_not_bp_skin))


indexM1 <- 1
indexM2 <- 2
log(topmodel_res[indexM1,pvars]/topmodel_res[indexM2,pvars])

selectorM1 <- c(1, topmodel_res[indexM1,(1:pvars-1)])
selectorM2 <- c(1, topmodel_res[indexM2,(1:pvars-1)])

source("~/R_programming/distbayesianmc/params_simulation/params_logit.R")
stan_code <- readChar(fileName, file.info(fileName)$size)
mod <- stan_model(model_code = stan_code, auto_write = T)

dataset_loaded <- f_dataset_loader("pima", nobs = nobs)

splitted_dataM1 <- f_pack_split_data(dataset_loaded$X[,selectorM1==1], dataset_loaded$y, ssplits=1, iseed=1, typesplit="random")
splitted_dataM1 <- f_prep_prior_logistic(splitted_dataM1, scale = 1)
res_approx_full <-  f_stan_sampling_splitted_data(mod, splitted_dataM1, dataset = "pima", i_seed = 1, iter = 1, typesplit = "random", nchain = 10000, typeprior="normal")
logbf_full <- res_approx_full$normconstcombined


splitted_dataM2 <- f_pack_split_data(dataset_loaded$X[,selectorM2==1], dataset_loaded$y, ssplits=1, iseed=1, typesplit="random")
splitted_dataM2 <- f_prep_prior_logistic(splitted_dataM2, scale = 1)
res_approx_small <-  f_stan_sampling_splitted_data(mod, splitted_dataM2, dataset = "pima", i_seed = 1, iter = 1, typesplit = "random", nchain = 10000, typeprior="normal")
logbf_small <-  res_approx_small$normconstcombined

logbf_small - logbf_full

      n.mil=1,
      seed=1,
      extra.arguments = list("Adaption_Iterations" = 2e5, 
                             "AlphaPriorSd" = 1,
                             "GaussianResidualPriorFamily" = 1,
                             "GaussianResidualPrior_UnifArg1" = 0.99,
                             "GaussianResidualPrior_UnifArg2" = 1.
                             )
    )
TopModels(pima.results)
ManhattanPlot(pima.results)
