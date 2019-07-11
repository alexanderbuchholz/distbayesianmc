// data {
//   int<lower=1> N;
//   real x[N];
// }
// transformed data {
//   matrix[N, N] K = cov_exp_quad(x, 1.0, 1.0);
//   vector[N] mu = rep_vector(0, N);
//   for (n in 1:N)
//     K[n, n] = K[n, n] + 0.1;
// }
// parameters {
//   vector[N] y;
// }
// model {
//   y ~ multi_normal(mu, K);
//}

data {
  int<lower=1> N;
  int<lower=0> P;
  real x[N];
  vector[N] y;
  real<lower=0> sigma; // the hyperparameter that gets adjusted
}
transformed data {
  vector[N] mu = rep_vector(0, N);
}
parameters {
  real beta[3];
  //real<lower=0> rho;
  //real<lower=0> alpha;
  //real<lower=0> sigma;
}
model {
  matrix[N, N] L_K;
  real sq_rho = square(beta[1]);
  matrix[N, N] K = cov_exp_quad(x, beta[2], sq_rho);
  real sq_sigma = square(beta[3]);

  // diagonal elements
  for (n in 1:N)
    K[n, n] = K[n, n] + sq_sigma;

  L_K = cholesky_decompose(K);

  //rho ~ inv_gamma(5, 5);
  beta[1] ~ normal(0, sigma^0.5);
  beta[2] ~ normal(0, sigma^0.5);
  beta[3] ~ normal(0, sigma^0.5);

  y ~ multi_normal_cholesky(mu, L_K);
}