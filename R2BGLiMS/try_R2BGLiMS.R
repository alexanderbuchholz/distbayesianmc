# R2BGLiMS
library(R2BGLiMS)
library(distbayesianmc)

dataset_loaded <- f_dataset_loader("pima")
pvars <- dim(dataset_loaded$X)[2]
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
