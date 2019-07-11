data {
  // Define variables in data
  // Number of observations (an integer)
  int<lower=0> N;
  // Number of parameters
  int<lower=0> P;
  // Variables
  int y[N];
  row_vector[P] x[N];
  real sigma; // the hyperparameter that gets adjusted
}

parameters {
  // Define parameters to estimate
  //vector[P] beta;
  vector[P] beta;
  //real alpha[1];
}


model {
  for (p in 1:P) {
    //beta[p] ~ normal(0, sigma^0.5);
    //target += normal_lpdf(beta[p] | 0, sigma^0.5);
    target += double_exponential_lpdf(beta[p] | 0, sigma);
  }
  for (n in 1:N){
    //target += bernoulli_lpmf( y[n]| inv_logit(x[n]*beta));
    target += bernoulli_logit_lpmf( y[n]| x[n]*beta);
    //y[n] ~ bernoulli(inv_logit(x[n]*beta));
  }
  
  
}