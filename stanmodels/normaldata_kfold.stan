data {
  int<lower=0> N;
  vector[N] x;
  
  int<lower = 0, upper = 1> loo; 
  vector[loo ? 3 : 0] x_loo;
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
  vector[3] loo_log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(x[n] | mu, sigma);
  if (loo) {
    for (n in 1:3) {
      loo_log_lik[n] = normal_lpdf(x_loo[n] | mu,sigma);
    }
  }
  else {
    for (n in 1:3) {
      loo_log_lik[n] = 0;
    }
  }
}
