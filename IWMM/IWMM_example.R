
source('IWMM/IWMM.R')


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Very simple example model
# Gaussian 1-dimensional data with unknown location and scale

stancode <- "data {
  int<lower=0> N;
  vector[N] x;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  x ~ normal(mu,sigma);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(x[n] | mu, sigma);
}"

stanmodel <- stan_model(model_code = stancode)




# Let us generate data from the true data generating mechanism
# Except set 1 observation by hand as an outlier



SEED <- 24
set.seed(SEED)
data_var = 1
n = as.integer(30)
x = rnorm(n = n, mean = 0,sd = data_var)
x[1] <- 15

# fit a model with all the observations
# We will use this model's posterior as the basis for importance sampling
standata = list(N = n, x = x)
fit <- sampling(stanmodel, data = standata, chains = 4, iter = 2000)

# Also fit a model without the outlier
# and compute the posterior mean
fit_loo <- sampling(stanmodel, data = list(N = n - 1, x = x[-c(1)]), chains = 4, iter = 2000)
(postmean_loo <- colMeans(as.matrix(fit_loo)[,1:2]))

# set the obs_weights for IWMM
# we are interested in the model without the outlier
# so we set the first weight to zero
obs_weights <- c(0,rep(1,n - 1))


# first without moment matching
iw <- IW(fit, obs_weights, expectation_fun =  function(x, upars, ...) {
  pars <- constrain_pars_stanfit(x, upars)
  pars[,1:2]},
  log_lik_fun = log_lik_stanfit, post_draws_fun = post_draws_stanfit,
  unconstrain_pars_fun = unconstrain_pars_stanfit)

(postmean_loo_is <- iw$expectation)

# the result is very different from the true mean





iw2 <- IW(fit, obs_weights, expectation_fun =  function(x, upars, ...) {
  pars <- constrain_pars_stanfit(x, upars)
  pars[,1:2]},
  log_lik_fun = log_lik_stanfit, post_draws_fun = post_draws_stanfit,
  unconstrain_pars_fun = unconstrain_pars_stanfit,
  moment_match = TRUE, k_threshold = 0.5,
  log_prob_upars_fun = log_prob_upars_stanfit,
  log_lik_upars_fun = log_lik_upars_stanfit)

(postmean_loo_is2 <- iw2$expectation)

# with moment matching the posterior mean is much closer to postmean_loo




# let us try leave-one-out cross-validation
# define the expectation function as the log likelihood
# of the left-out observation
expectation_fun_elpd_log <- function(x, upars, i, ...) {
  matrix(log_lik_upars_stanfit(x, upars, ...)[,i])
}



loo_iw <- IW(x = fit, obs_weights, expectation_fun = expectation_fun_elpd_log,
                  log_lik_fun = log_lik_stanfit, post_draws_fun = post_draws_stanfit, log_expectation_fun = TRUE, i = 1,
                  moment_match = TRUE, k_threshold = 0.5,
                  unconstrain_pars_fun = unconstrain_pars_stanfit,
                  log_prob_upars_fun = log_prob_upars_stanfit,
                  log_lik_upars_fun = log_lik_upars_stanfit,
                  split = TRUE, restart_transform = TRUE)


(loo_is <- log(loo_iw$expectation))


# let us compare the result to the loo package
looobj = rstan::loo(fit, moment_match = TRUE, k_threshold = 0.5)
(looobj$pointwise[1])

# the result is similar



