data {
  int<lower=0> N;
  int<lower=0> J;
  vector[N] y;
  vector[N] x;
  int<lower=1,upper=J> county[N];
}
parameters {
  real<lower=0> sigma_y;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  vector[J] a;
  vector[J] b;
  real mu_a;
  real mu_b;
}
transformed parameters {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] = a[county[i]] + x[i] * b[county[i]];
}
model {
  sigma_y ~ normal(0, 1);
  sigma_a ~ normal(0, 1);
  sigma_b ~ normal(0, 1);

  mu_a ~ normal(0, 10);
  mu_b ~ normal(0, 10);

  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ normal(y_hat, sigma_y);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | y_hat[n] , sigma_y);
}

