data {
  int<lower=1> K;
  int<lower=1> N;
  matrix[N,K] x;
  int y[N];
  vector[N] offset;

  real beta_prior_scale;
  real alpha_prior_scale;

  // whether to predict leave-one-out obs. or not
  int<lower = 0, upper = 1> loo; 
  vector[loo ? K : 0] x_loo;
  int y_loo[loo ? 1 : 0]; 
  real offset_loo[loo ? 1 : 0]; 
}
parameters {
  vector[K] beta;
  real intercept;
}
model {
  y ~ poisson(exp(x * beta + intercept + offset));
  for (k in 1:K) {
    beta[k] ~ normal(0,beta_prior_scale);
  }
  intercept ~ normal(0,alpha_prior_scale);
}
generated quantities {
  vector[N] log_lik;
  real loo_log_lik; 
  for (n in 1:N)
    log_lik[n] = poisson_lpmf(y[n] | exp(x[n] * beta + intercept + offset[n]));
  if (loo) {
    loo_log_lik = poisson_lpmf(y_loo | exp(dot_product(x_loo , beta) + intercept + offset_loo[1]));
  }
  else {
    loo_log_lik = 0;
  }

}
