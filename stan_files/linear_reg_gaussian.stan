data {
  int<lower=1> N; // Number of data
  int<lower=1> P; // Number of covariates
  matrix[N, P] x;
  real y[N];
  real<lower=0> sigma; // the hyperparameter that gets adjusted
}

parameters {
  vector[P+2] beta;
  //real alpha;
  //real<lower=0> sigma_x;
}

model {
  // Strongly regularizing priors
  //beta[1:P] ~ double_exponential(0, sigma);
  
  for (p in 1:P) {
    //beta[p] ~ normal(0, sigma^0.5);
    //target += normal_lpdf(beta[p] | 0, sigma^0.5);
    target += normal_lpdf(beta[p] | 0, sigma^0.5);
  }
  //beta[1:P] ~ normal(0, sigma);
  //beta[P+1] ~ normal(0, sigma^0.5);
  target += normal_lpdf(beta[P+1] | 0, sigma^0.5);
  target += normal_lpdf(beta[P+2] | 0, sigma^0.5);
  //beta[P+2] ~ normal(0, sigma);
  //alpha ~ normal(0, 2);
  //sigma ~ normal(0, 2);
  
  //y ~ normal((x * beta[1:P]) + beta[P+1], exp(beta[P+2]));
  //y ~ normal((x * beta[1:P]) + beta[P+1], 1);
  target += normal_lpdf(y | (x * beta[1:P]) + beta[P+1], exp(beta[P+2]));
}
