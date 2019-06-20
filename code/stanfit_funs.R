# create a named list of draws for use with rstan methods
.rstan_relist <- function (x, skeleton) {
  out <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) {
    dim(out[[i]]) <- dim(skeleton[[i]])
  }
  out
}

# rstan helper function to get dims of parameters right
.create_skeleton <- function (pars, dims) {
  out <- lapply(seq_along(pars), function(i) {
    len_dims <- length(dims[[i]])
    if (len_dims < 1) return(0)
    return(array(0, dim = dims[[i]]))
  })
  names(out) <- pars
  out
}

# extract original posterior draws
post_draws_stanfit <- function(x, ...) {
  as.matrix(x)
}

# compute a vector of log-likelihood values for the ith observation
log_lik_stanfit <- function(x, i, parameter_name = "log_lik", ...) {
  loo::extract_log_lik(x, parameter_name, merge_chains = FALSE)[, , i]
}

# transform parameters to the unconstraint space
unconstrain_pars_stanfit <- function(x, pars, ...) {
  skeleton <- .create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
  upars <- apply(pars, 1, FUN = function(theta) {
    rstan::unconstrain_pars(x, .rstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(upars))) {
    dim(upars) <- c(1, length(upars))
  }
  t(upars)
}

# compute log_prob for each posterior draws on the unconstrained space
log_prob_upars_stanfit <- function(x, upars, ...) {
  apply(upars, 1, rstan::log_prob, object = x,
        adjust_transform = TRUE, gradient = FALSE)
}

# compute log_lik values based on the unconstrained parameters
log_lik_upars_stanfit <- function(x, upars, i, parameter_name = "log_lik",
                                  ...) {
  S <- nrow(upars)
  out <- numeric(S)
  for (s in seq_len(S)) {
    out[s] <- rstan::constrain_pars(x, upars = upars[s, ])[[parameter_name]][i]
  }
  out
}

#' Leave-one-out cross-validation
#' 
#' The \code{loo} method for stanfit objects ---a wrapper around the
#' \code{loo.array} method from the \pkg{loo} package--- computes approximate
#' leave-one-out cross-validation using Pareto smoothed importance sampling
#' (PSIS-LOO CV).
#' 
#' @param x stanfit object
#' @param pars Name of parameter, transformed parameter, or generated quantity in
#'   the Stan program corresponding to the pointwise log-likelihood. If not
#'   specified the default behavior is to look for \code{"log_lik"}.
#' @param ... Further arguments passed to underlying methods.
#' @param cores Number of cores to use for parallelization. Passed to
#'   \code{\link[loo]{loo.array}}.
#' @param save_psis Should the intermediate results from \code{\link[loo]{psis}}
#'   be saved in the returned object? The default is \code{FALSE}. This can be
#'   useful to avoid repeated computation when using other functions in the
#'   \pkg{loo} and \pkg{bayesplot} packages.
#' @param k_thres Threshold value for Pareto k values above which
#'   the moment matching algorithm is used.
#' 
#' @details Stan does not automatically compute and store the log-likelihood. It
#'   is up to the user to incorporate it into the Stan program if it is to be
#'   extracted after fitting the model. In a Stan model, the pointwise log
#'   likelihood can be coded as a vector in the transformed parameters block
#'   (and then summed up in the model block) or it can be coded entirely in the
#'   generated quantities block. We recommend using the generated quantities
#'   block so that the computations are carried out only once per iteration
#'   rather than once per HMC leapfrog step.
#'
#'   For example, the following is the \code{generated quantities} block for
#'   computing and saving the log-likelihood for a linear regression model with
#'   \code{N} data points, outcome \code{y}, predictor matrix \code{X},
#'   coefficients \code{beta}, and standard deviation \code{sigma}:
#'
#' \code{vector[N] log_lik;}
#'
#' \code{for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | X[n, ] * beta, sigma);}
#'
loo_stanfit <-
  function(x,
           pars = "log_lik",
           ...,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1),
           k_thres = 0.5) {
    stopifnot(length(pars) == 1L)
    LLarray <- loo::extract_log_lik(stanfit = x,
                                    parameter_name = pars,
                                    merge_chains = FALSE)
    r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
    
    loo <- SW(loo::loo.array(LLarray,
                          r_eff = r_eff,
                          cores = cores,
                          save_psis = save_psis))
    
    # if there are high Pareto k values, use moment matching algorithm
    if (any(loo$diagnostics$pareto_k >= k_thres)) {
      loo <- mm_loo_stanfit(
        x, loo = loo, k_thres = k_thres, 
        cores = cores, parameter_name = pars, ...
      )
    }
    # Warn if some Pareto ks are still high
    psislw_warnings(loo$diagnostics$pareto_k)
    
    loo
  }

# wrapper around mm_loo_manual
mm_loo_stanfit <- function(x, loo, ...) {
  mm_loo_manual(
    x, loo = loo, 
    post_draws = post_draws_stanfit, 
    log_lik = log_lik_stanfit, 
    unconstrain_pars = unconstrain_pars_stanfit,
    log_prob_upars = log_prob_upars_stanfit,
    log_lik_upars = log_lik_upars_stanfit,
    ...
  )
}
