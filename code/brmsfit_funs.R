# extract original posterior draws
post_draws_brmsfit <- function(x, ...) {
  as.matrix(x, ...)
}

# compute a vector of log-likelihood values for the ith observation
log_lik_brmsfit <- function(x, i, ...) {
  as.vector(log_lik(x, newdata = x$data[i, , drop = FALSE], ...))
}

# transform parameters to the unconstraint space
unconstrain_pars_brmsfit <- function(x, pars, ...) {
  unconstrain_pars_stanfit(x$fit, pars = pars, ...)
}

# compute log_prob for each posterior draws on the unconstrained space
log_prob_upars_brmsfit <- function(x, upars, ...) {
  log_prob_upars_stanfit(x$fit, upars = upars, ...)
}

# transform parameters to the constraint space
update_pars_brmsfit <- function(x, upars, ...) {
  # list with one element per posterior draw
  pars <- apply(upars, 1, rstan::constrain_pars, object = x$fit)
  # transform samples
  nsamples <- length(pars)
  pars <- unlist(pars)
  npars <- length(pars) / nsamples
  dim(pars) <- c(npars, nsamples)
  # add dummy 'lp__' samples
  pars <- rbind(pars, rep(0, nsamples))
  # bring samples into the right structure
  new_samples <- brms:::named_list(
    x$fit@sim$fnames_oi_old, list(numeric(nsamples))
  )
  stopifnot(length(new_samples) == nrow(pars))
  for (i in seq_len(npars)) {
    new_samples[[i]] <- pars[i, ]
  }
  # create new sim object to overwrite x$fit@sim
  x$fit@sim <- list(
    samples = list(new_samples),
    iter = nsamples,
    thin = 1,
    warmup = 0,
    chains = 1,
    n_save = nsamples,
    warmup2 = 0,
    permutation = list(seq_len(nsamples)),
    pars_oi = x$fit@sim$pars_oi_old,
    dims_oi = x$fit@sim$dims_oi_old,
    fnames_oi = x$fit@sim$fnames_oi_old,
    n_flatnames = length(x$fit@sim$fnames_oi_old)
  ) 
  x$fit@stan_args <- list(
    list(chain_id = 1, iter = nsamples, thin = 1, warmup = 0)
  )
  brms:::rename_pars(x)
}

# compute log_lik values based on the unconstrained parameters
log_lik_upars_brmsfit <- function(x, upars, i, samples = NULL, 
                                  subset = NULL, ...) {
  # do not pass subset or nsamples further to avoid subsetting twice
  x <- update_pars_brmsfit(x, upars = upars, ...)
  log_lik_brmsfit(x, i = i, ...)
}

# wrapper around mm_loo_manual
mm_loo_brmsfit <- function(x, loo, ...) {
  mm_loo_manual(
    x, loo = loo, 
    post_draws = post_draws_brmsfit, 
    log_lik = log_lik_brmsfit, 
    unconstrain_pars = unconstrain_pars_brmsfit,
    log_prob_upars = log_prob_upars_brmsfit,
    log_lik_upars = log_lik_upars_brmsfit, 
    ...
  )
}

