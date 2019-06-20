library(rstan)
library(brms)
source("rstan_code/mm_loo_funs.R")
source("rstan_code/stanfit_funs.R")
source("rstan_code/brmsfit_funs.R")
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

CHAINS <- 4
ITER <- 2000
SEED <- 12
set.seed(SEED)

# generate toy data
data_sd <- 1
n <- as.integer(30)
x <- rnorm(n=n, mean=0, sd=data_sd)
x_tilde <- 11
x[1] <- x_tilde

# fit the model with rstan
standata <- list(N = n, x = x, loo = 0, x_loo = numeric(0))
fit_rstan <- stan(file = 'stanmodels/normaldata.stan', data = standata,
            chains = CHAINS, seed = SEED, iter = ITER,
            control = list(adapt_delta = 0.95), save_warmup = FALSE)
print(fit_rstan)

loo_rstan <- loo(fit_rstan,save_psis = TRUE)
plot(loo_rstan)

loo_rstan2 <- mm_loo_stanfit(fit_rstan, loo = loo_rstan, 
                             split = TRUE, save_psis = TRUE)

loo_rstan3 <- loo_stanfit(fit_rstan, split = TRUE, save_psis = TRUE)

loo_rstan$estimates
loo_rstan2$estimates
loo_rstan3$estimates


# fit the model with brms
fit_brms <- brm(x ~ 1, data = data.frame(x), save_all_pars = TRUE,
                chains = CHAINS, seed = SEED, iter = ITER,
                control = list(adapt_delta = 0.95))
print(fit_brms)

loo_brms <- loo(fit_brms)
plot(loo_brms)

loo_brms2 <- mm_loo_brmsfit(fit_brms, loo = loo_brms, split = TRUE)

loo_brms$estimates
loo_brms2$estimates



