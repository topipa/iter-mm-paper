
library(rstan)
rstan_options(auto_write = TRUE)
cores <- 4
options(mc.cores = cores)
source('IWMM.R')

stancode <- "data {
  int<lower=0> N;
  vector[N] x;
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
  for (n in 1:N)
    log_lik[n] = normal_lpdf(x[n] | mu, sigma);
}"


stanmodel <- stan_model(model_code = stancode)


SEED <- 21
set.seed(SEED)

# generate toy data
data_var = 1 # std
n = as.integer(200)
x = rnorm(n = n, mean = 0,sd = data_var)

standata = list(N = n, x = x)
fit <- sampling(stanmodel, data = standata, chains = 4, iter = 4000)

samples <- as.data.frame(fit)
S <- nrow(samples)
post <- samples[,1:2]
ll <- samples[,4:(4 + n - 1)]


bb_S <- 500


bb_postmeans_is <- matrix(0,bb_S,2)
ks <- rep(0,bb_S)

bb_postmeans_is_mm <- matrix(0,bb_S,2)
ks_mm <- rep(0,bb_S)



iw_list <- parallel::mclapply(X = seq(bb_S), mc.cores = cores, FUN = function(i) {
  obs_w <- c(MCMCpack::rdirichlet(1,rep(1,n)))
  obs_w <- obs_w/sum(obs_w)*n

  obs_w_matrix <- matrix(obs_w,S,n, byrow = TRUE)
  lw <- rowSums((obs_w_matrix - 1) * ll)
  lw <- lw - matrixStats::logSumExp(lw)
  w <- exp(lw)
  psis <- loo::psis(lw)



  iw <- IW(fit, obs_w, expectation_fun =  function(x, upars, ...) {
    # omitted constraining because model does not have constrained params
    pars <- upars
    pars[,1:2]},
    log_lik_fun = log_lik_stanfit, post_draws_fun = post_draws_stanfit,
    unconstrain_pars_fun = unconstrain_pars_stanfit,
    moment_match = TRUE, k_threshold = 0.5,
    log_prob_upars_fun = log_prob_upars_stanfit,
    log_lik_upars_fun = log_lik_upars_stanfit)

  c(colSums(w * post),psis$diagnostics$pareto_k,iw$expectation,iw$pareto_k)
})

iw_list <- array(as.numeric(unlist(iw_list)),c(6,bb_S))

bb_postmeans_is <- t(iw_list[1:2,])
ks <- iw_list[3,]
bb_postmeans_is_mm <- t(iw_list[4:5,])
ks_mm <- iw_list[6,]





plot(post$mu,post$logsigma)
points(bb_postmeans_is[,1], bb_postmeans_is[,2], col = 'red')
points(bb_postmeans_is_mm[,1], bb_postmeans_is_mm[,2], col = 'blue')


orders <- order(ks)
plot(ks[orders])
points(ks_mm[orders], col = 'blue')

matrixStats::colVars(as.matrix(post))
matrixStats::colVars(bb_postmeans_is)
matrixStats::colVars(bb_postmeans_is_mm)



