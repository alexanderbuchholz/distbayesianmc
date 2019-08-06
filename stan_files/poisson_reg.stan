// negative binomial regression

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
    target += normal_lpdf(beta[p] | 0, sigma^0.5);
  }
  for (n in 1:N){
    target += poisson_lpmf( y[n]| Phi_approx(x[n]*beta));
    //target += bernoulli_logit_lpmf( y[n]| x[n]*beta);
  }
  
  
}