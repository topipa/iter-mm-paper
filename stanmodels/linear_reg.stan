data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] x;
  vector[N] y;
  real<lower=0> beta_prior_scale;
  real<lower=0> alpha_prior_scale;
  
  // whether to predict leave-one-out obs. or not
  int<lower = 0, upper = 1> loo; 
  vector[loo ? K : 0] x_loo;
  real y_loo[loo ? 1 : 0]; 

}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma; // error scale
}
model {
  y ~ normal(x * beta + alpha, sigma); // likelihood
  beta ~ normal(0,beta_prior_scale);
  alpha ~ normal(0,alpha_prior_scale);
  sigma ~ exponential(1);
}
generated quantities {
  vector[N] log_lik;
  real loo_log_lik; 
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | x[n] * beta + alpha, sigma);
  if (loo) {
    loo_log_lik = normal_lpdf(y_loo | dot_product(x_loo , beta) + alpha, sigma);
  }
  else {
    loo_log_lik = 0;
  }
}

