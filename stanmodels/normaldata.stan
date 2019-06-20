data {
  int<lower=0> N;
  vector[N] x;
  
  int<lower = 0, upper = 1> loo; 
  vector[loo ? 1 : 0] x_loo;
}
parameters {
  real mu;
  real logsigma;
}
transformed parameters {
  real<lower=0> sigma;
  sigma = exp(logsigma);
}
model {
  x ~ normal(mu,sigma);
}
generated quantities {
  vector[N] log_lik;
  real loo_log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(x[n] | mu, sigma);
  if (loo) {
    loo_log_lik = normal_lpdf(x_loo | mu,sigma);
  }
  else {
    loo_log_lik = 0;
  }
}
