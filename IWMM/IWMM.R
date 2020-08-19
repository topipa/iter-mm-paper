#' Generic importance weighted moment matching (IWMM) importance sampling algorithm.
#' for cross-validation, bootstrap or similar purposes. Can compute arbitrary expectations
#' over arbitrary distributions that are variations of the model posterior distribution.
#'
#'
#' @param x A fitted model object.
#' @param obs_weights A vector with length equal to number of observations in
#' the model x. The weights represent the target density for importance sampling.
#' Each observation is given a nonnegative weight. A vector of all ones
#' represents the original posterior distribution. A weight zero represents omitting
#' a certain observation.
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `x` and `upars`. Must return a matrix with
#' rows equal to the number of draws.
#' @param log_lik_fun A function that takes `x` and returns a matrix of log-likelihood draws
#'  based on the model `x`. The rows indicate draws and columns indicate observations.
#' @param post_draws_fun A function that takes `x` as the first argument and returns
#' a matrix of posterior draws of the model parameters.
#' @param unconstrain_pars_fun A function that takes arguments `x`, and `pars` and
#' returns posterior draws on the unconstrained space based on the posterior
#' draws on the constrained space passed via `pars`.
#' @param moment_match Logical indicating whether to use moment matching
#' to improve the accuracy of importance sampling. Defaults to FALSE.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @param ... Further arguments passed to the custom functions documented above.
#'
#' @return Returns a list with 3 elements: the expectation of interest,
#' and 2 pareto k diagnostic values.
IW <- function(x, obs_weights, expectation_fun, log_lik_fun, post_draws_fun, unconstrain_pars_fun, moment_match = FALSE, log_expectation_fun = FALSE, ...) {

  pars <- post_draws_fun(x, ...)
  # transform the model parameters to unconstrained space
  upars <- unconstrain_pars_fun(x, pars = pars, ...)
  log_lik <- log_lik_fun(x, ...)
  S <- nrow(log_lik)
  N <- ncol(log_lik)

  obs_weights <- c(obs_weights)
  obs_weights <- matrix(obs_weights, S, N, byrow = TRUE)

  lw <- rowSums((obs_weights - 1) * log_lik)
  lw_psis <- suppressWarnings(loo::psis(lw))
  lw <- as.vector(weights(lw_psis))
  k <- lw_psis$diagnostics$pareto_k


  if (moment_match) {
    IWMM(x, obs_weights, upars, lw_psis, expectation_fun, log_expectation_fun, ...)
  }
  else {
    if (log_expectation_fun) {
      expectation <- exp(matrixStats::colLogSumExps(lw + expectation_fun(x, upars, ...)))
    }
    else {
      w <- exp(lw)
      expectation <- colSums(w * expectation_fun(x, upars, ...))
    }
    lwf_psis <- suppressWarnings(loo::psis(expectation))
    kf <- lwf_psis$diagnostics$pareto_k


    list("expectation" = expectation, "pareto_k" = k, "pareto_kf" = kf)
  }


}




# restart_transform ignored if split is FALSE
#' @param x A fitted model object.
#' @param obs_weights A vector with length equal to number of observations in
#' the model x. The weights represent the target density for importance sampling.
#' Each observation is given a nonnegative weight. A vector of all ones
#' represents the original posterior distribution. A weight zero represents omitting
#' a certain observation.
#' @param upars A matrix of model parameters in unconstrained space.
#' @param lw_psis A psis object from the loo package computed from importance weights
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `x` and `upars`.
#' @param log_prob_upars_fun A function that takes `x` and `upars` and returns a
#' vector of model posterior density values based on `upars`.
#' @param log_lik_upars_fun A function that takes `x` and `upars` and returns a matrix of log-likelihood draws
#' based on the model `x` and `upars`. The rows indicate draws and columns indicate observations.
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @param k_threshold Threshold for pareto k values above which the moment matching is used.
#' Set to 0.5 by default. Lower (higher) threshold leads
#' to better (worse) accuracy but slower (faster) computation.
#' @param cov_transform Logical; Indicate whether to match the covariance matrix of the
#'   samples or not. If `FALSE`, only the mean and marginal variances are
#'   matched.
#' @param split Logical; Indicate whether to do the split transformation or not
#' at the end of moment matching. FALSE by default.
#' @param restart_transform Logical; When split is TRUE, indicates whether to
#' start the second transformation from the original model parameters
#' or the transformed parameters.
#' @param ... Further arguments passed to the custom functions documented above.
#'
#' @return Returns a list with 3 elements: the expectation of interest,
#' and 2 pareto k diagnostic values.
IWMM <- function(x, obs_weights, upars, lw_psis, expectation_fun, log_prob_upars_fun, log_lik_upars_fun, log_expectation_fun, k_threshold = 0.5, cov_transform = FALSE, split = FALSE, restart_transform = FALSE, ...) {
  upars_orig <- upars
  orig_log_prob <- log_prob_upars_fun(x, upars = upars, ...)

  npars <- ncol(upars)
  S <- nrow(upars)
  cov_transform <- cov_transform && S >= 10 * npars

  lw <- as.vector(weights(lw_psis))
  k <- lw_psis$diagnostics$pareto_k
  lw_orig <- lw

  # initialize objects that keep track of the total transformation
  total_shift <- rep(0, npars)
  total_scaling <- rep(1, npars)
  total_mapping <- diag(npars)


  while (k > k_threshold) {


    # 1. match means
    trans <- shift(x, upars, lw)
    quantities <- update_quantities(
      x, obs_weights,
      upars = trans$upars,
      orig_log_prob = orig_log_prob,
      log_prob_upars = log_prob_upars_fun,
      log_lik_upars = log_lik_upars_fun,
      ...
    )
    if (quantities$k < k) {
      upars <- trans$upars
      total_shift <- total_shift + trans$shift

      lw <- quantities$lw
      k <- quantities$k
      next
    }

    # 2. match means and marginal variances
    trans <- shift_and_scale(x, upars, lw)
    quantities <- update_quantities(
      x, obs_weights,
      upars = trans$upars,
      orig_log_prob = orig_log_prob,
      log_prob_upars = log_prob_upars_fun,
      log_lik_upars = log_lik_upars_fun,
      ...
    )
    if (quantities$k < k) {
      upars <- trans$upars
      total_shift <- total_shift + trans$shift
      total_scaling <- total_scaling * trans$scaling

      lw <- quantities$lw
      k <- quantities$k
      next
    }

    # 3. match means and covariances
    trans <- shift_and_cov(x, upars, lw)
    quantities <- update_quantities(
      x, obs_weights,
      upars = trans$upars,
      orig_log_prob = orig_log_prob,
      log_prob_upars = log_prob_upars_fun,
      log_lik_upars = log_lik_upars_fun,
      ...
    )
    if (quantities$k < k) {
      upars <- trans$upars
      total_shift <- total_shift + trans$shift
      total_mapping <- trans$mapping %*% total_mapping

      lw <- quantities$lw
      k <- quantities$k
      next
    }


    break
  }

  # prepare for split and check kfs

  if (restart_transform) {
    upars2 <- upars_orig
    total_shift2 <- rep(0, npars)
    total_scaling2 <- rep(1, npars)
    total_mapping2 <- diag(npars)
    lw <- lw_orig

  }
  else {
    upars2 <- upars
    # initialize objects that keep track of the total transformation
    total_shift2 <- total_shift
    total_scaling2 <- total_scaling
    total_mapping2 <- total_mapping
  }

  if (log_expectation_fun) {
    lwf <- lw + expectation_fun(x, upars2, ...)
    # Here we expect that expectation_fun is nonnegative (before log)')
  }
  else {
    lwf <- lw + log(abs(expectation_fun(x, upars2, ...)))
  }
  psisf <- suppressWarnings(loo::psis(lwf))
  kf <- psisf$diagnostics$pareto_k


  if (split) {


    if (ncol(lwf) > 1) {
      stop('Using split = TRUE is not yet supported for expectation functions that return a matrix.
           As a workaround, you can wrap your function call using apply.')
    }

    lwf <- as.vector(weights(psisf))

    while (kf > k_threshold) {

      # 1. match means
      trans <- shift(x, upars2, lwf)
      quantities <- update_quantities2(
        x, obs_weights,
        upars2 = trans$upars,
        orig_log_prob = orig_log_prob,
        log_prob_upars = log_prob_upars_fun,
        log_lik_upars = log_lik_upars_fun,
        expectation_fun,
        log_expectation_fun = log_expectation_fun,
        ...
      )
      if (quantities$kf < kf) {
        upars2 <- trans$upars
        total_shift2 <- total_shift2 + trans$shift

        lwf <- quantities$lwf
        kf <- quantities$kf
        next
      }

      # 2. match means and variances
      trans <- shift_and_scale(x, upars2, lwf)
      quantities <- update_quantities2(
        x, obs_weights,
        upars2 = trans$upars,
        orig_log_prob = orig_log_prob,
        log_prob_upars = log_prob_upars_fun,
        log_lik_upars = log_lik_upars_fun,
        expectation_fun,
        log_expectation_fun = log_expectation_fun,
        ...
      )
      if (quantities$kf < kf) {
        upars2 <- trans$upars
        total_shift2 <- total_shift2 + trans$shift
        total_scaling2 <- total_scaling2 * trans$scaling

        lwf <- quantities$lwf
        kf <- quantities$kf
        next
      }

      # 3. match means and covariances
      trans <- shift_and_cov(x, upars2, lwf)
      quantities <- update_quantities2(
        x, obs_weights,
        upars2 = trans$upars,
        orig_log_prob = orig_log_prob,
        log_prob_upars = log_prob_upars_fun,
        log_lik_upars = log_lik_upars_fun,
        expectation_fun,
        log_expectation_fun = log_expectation_fun,
        ...
      )
      if (quantities$kf < kf) {
        upars2 <- trans$upars
        total_shift2 <- total_shift2 + trans$shift
        total_mapping2 <- trans$mapping %*% total_mapping2

        lwf <- quantities$lwf
        kf <- quantities$kf
        next
      }

      break
    }

    # second trasnsformations are done
    # we have updated upars2, lwf and kf
    # now compute split weights and split pars

    S <- nrow(upars2)
    S_half <- as.integer(0.5 * S)


    # original parameters
    mean_original <- colMeans(upars_orig)

    # accumulated affine transformation
    upars_T1 <- sweep(upars_orig, 2, mean_original, "-")
    upars_T1 <- sweep(upars_T1, 2, total_scaling, "*")
    if (cov_transform) {
      upars_T1 <- tcrossprod(upars_T1, total_mapping)
    }
    upars_T1 <- sweep(upars_T1, 2, total_shift + mean_original, "+")

    upars_T2 <- sweep(upars_orig, 2, mean_original, "-")
    upars_T2 <- sweep(upars_T2, 2, total_scaling2, "*")
    if (cov_transform) {
      upars_T2 <- tcrossprod(upars_T2, total_mapping2)
    }
    upars_T2 <- sweep(upars_T2, 2, total_shift2 + mean_original, "+")



    # inverse accumulated affine transformation
    upars_T2_T1inv <- sweep(upars_T2, 2, mean_original + total_shift2, "-")
    if (cov_transform) {
      upars_T2_T1inv <- tcrossprod(upars_T2_T1inv, solve(total_mapping))
    }
    upars_T2_T1inv <- sweep(upars_T2_T1inv, 2, total_scaling, "/")
    upars_T2_T1inv <- sweep(upars_T2_T1inv, 2, mean_original + total_shift2 - total_shift, "+")

    upars_T1_T2inv <- sweep(upars_T1, 2, mean_original + total_shift, "-")
    if (cov_transform) {
      upars_T1_T2inv <- tcrossprod(upars_T1_T2inv, solve(total_mapping2))
    }
    upars_T1_T2inv <- sweep(upars_T1_T2inv, 2, total_scaling2, "/")
    upars_T1_T2inv <- sweep(upars_T1_T2inv, 2, mean_original + total_shift - total_shift2, "+")




    # these are the real used draws
    # first half of upars_trans are T1(theta)
    # second half are T2(theta)
    upars_trans <- upars_T2
    take <- seq_len(S_half)
    upars_trans[take, ] <- upars_T1[take, , drop = FALSE]


    # then we need two sets of pseudo draws
    upars_trans_inv1 <- upars_T2_T1inv
    take <- seq_len(S)[-seq_len(S_half)]
    upars_trans_inv1[take, ] <- upars_orig[take, , drop = FALSE]

    upars_trans_inv2 <- upars_orig
    take <- seq_len(S)[-seq_len(S_half)]
    upars_trans_inv2[take, ] <- upars_T1_T2inv[take, , drop = FALSE]


    log_prob_trans <- log_prob_upars_fun(x, upars = upars_trans, ...)
    log_lik_trans <- log_lik_upars_fun(x, upars = upars_trans, ...)

    log_prob_trans_inv1 <- log_prob_upars_fun(x, upars = upars_trans_inv1, ...)
    log_prob_trans_inv2 <- log_prob_upars_fun(x, upars = upars_trans_inv2, ...)


    lw_trans <- rowSums((obs_weights - 1) * log_lik_trans) + log_prob_trans -
      log_prob_trans_inv1 - log(prod(total_scaling2))  - log(det(total_mapping2)) -
      log(1 + exp(log_prob_trans_inv2 - log(prod(total_scaling))  - log(det(total_mapping)) -
                    log_prob_trans_inv1 - log(prod(total_scaling2)) - log(det(total_mapping2))     ))


    lw_trans_psis <- suppressWarnings(loo::psis(lw_trans))
    lw_trans <- as.vector(weights(lw_trans_psis))


    # replace upars and lw
    lw <- lw_trans
    upars <- upars_trans



  }
  else {
    # if not splitting, warn about high pareto ks
    if (any(kf > k_threshold)) {
      warning('Importance sampling may be unreliable. Consider setting split to TRUE.')
    }
  }


  if (log_expectation_fun) {
    expectation <- exp(matrixStats::colLogSumExps(lw + expectation_fun(x, upars, ...)))
  }
  else {
    w <- exp(lw)
    expectation <- colSums(w * expectation_fun(x, upars, ...))
  }

  list("expectation" = expectation, "pareto_k" = k, "pareto_kf" = kf)
}






####################### HELPER FUNCTIONS

# extract original posterior draws
post_draws_stanfit <- function(x, ...) {
  as.matrix(x)
}

log_lik_stanfit <- function(x, i, parameter_name = "log_lik", ...) {
  loo::extract_log_lik(x, parameter_name, merge_chains = TRUE)
}

log_lik_i_stanfit <- function(x, i, parameter_name = "log_lik", ...) {
  loo::extract_log_lik(x, parameter_name, merge_chains = FALSE)[, , i]
}




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


constrain_pars_stanfit <- function(x, upars, ...) {
  npars <- rstan::get_num_upars(x)
  pars <- apply(upars, 1, FUN = function(theta){
    unlist(rstan::constrain_pars(x,theta))[seq(npars)]
  })
  # for one parameter models
  if (is.null(dim(pars))) {
    dim(pars) <- c(1, length(pars))
  }
  t(pars)
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


# create a named list of draws for use with rstan methods
.rstan_relist <- function (x, skeleton) {
  out <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) {
    dim(out[[i]]) <- dim(skeleton[[i]])
  }
  out
}

# compute log_prob for each posterior draws on the unconstrained space
log_prob_upars_stanfit <- function(x, upars, ...) {
  apply(upars, 1, rstan::log_prob, object = x,
        adjust_transform = TRUE, gradient = FALSE)
}



log_lik_upars_stanfit <- function(x, upars, parameter_name = "log_lik",
                                  ...) {
  temp <- rstan::constrain_pars(x, upars = upars[1, ])[[parameter_name]]
  n <- length(temp)
  S <- nrow(upars)
  out <- matrix(0,S,n)
  for (s in seq_len(S)) {
    out[s,] <- rstan::constrain_pars(x, upars = upars[s, ])[[parameter_name]]
  }
  out

}


shift <- function(x, upars, lw) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lw) * upars)
  shift <- mean_weighted - mean_original
  # transform posterior draws
  upars_new <- sweep(upars, 2, shift, "+")
  list(
    upars = upars_new,
    shift = shift
  )
}

shift_and_scale <- function(x, upars, lwi) {
  # compute moments using log weights
  S <- dim(upars)[1]
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  mii <- exp(lwi)* upars^2
  mii <- colSums(mii) - mean_weighted^2
  mii <- mii*S/(S-1)
  scaling <- sqrt(mii / matrixStats::colVars(upars))
  # transform posterior draws
  upars_new <- sweep(upars, 2, mean_original, "-")
  upars_new <- sweep(upars_new, 2, scaling, "*")
  upars_new <- sweep(upars_new, 2, mean_weighted, "+")

  list(
    upars = upars_new,
    shift = shift,
    scaling = scaling
  )
}

shift_and_cov <- function(x, upars, lwi, ...) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  covv <- stats::cov(upars)
  wcovv <- stats::cov.wt(upars, wt = exp(lwi))$cov
  chol1 <- tryCatch(
    {
      chol(wcovv)
    },
    error = function(cond)
    {
      return(NULL)
    }
  )
  if (is.null(chol1)) {
    mapping <- diag(length(mean_original))
  }
  else {
    chol2 <- chol(covv)
    mapping <- t(chol1) %*% solve(t(chol2))
  }
  # transform posterior draws
  upars_new <- sweep(upars, 2, mean_original, "-")
  upars_new <- tcrossprod(upars_new, mapping)
  upars_new <- sweep(upars_new, 2, mean_weighted, "+")
  colnames(upars_new) <- colnames(upars)

  list(
    upars = upars_new,
    shift = shift,
    mapping = mapping
  )
}


update_quantities <- function(x, obs_weights, upars, orig_log_prob,
                              log_prob_upars, log_lik_upars,
                              ...) {
  log_prob_new <- log_prob_upars(x, upars = upars, ...)
  log_lik_new <- log_lik_upars(x, upars = upars, ...)

  lw_new <- rowSums((obs_weights - 1) * log_lik_new) + log_prob_new - orig_log_prob

  psis_new <- suppressWarnings(loo::psis(lw_new))
  k_new <- psis_new$diagnostics$pareto_k
  lw_new <- as.vector(weights(psis_new))

  # gather results
  list(
    lw = lw_new,
    k = k_new
  )
}


update_quantities2 <- function(x, obs_weights, upars2, orig_log_prob,
                               log_prob_upars, log_lik_upars,
                               expectation_fun, log_expectation_fun,
                               ...) {


  log_prob_new <- log_prob_upars(x, upars = upars2, ...)
  log_lik_new <- log_lik_upars(x, upars = upars2, ...)
  lw_new <- rowSums((obs_weights - 1) * log_lik_new) + log_prob_new - orig_log_prob

  if (log_expectation_fun) {
    lwf_new <- lw_new + expectation_fun(x, upars2, ...)
  }
  else {
    lwf_new <- lw_new + log(abs(expectation_fun(x, upars2, ...)))
  }

  psisf_new <- suppressWarnings(loo::psis(lwf_new))
  kf_new <- psisf_new$diagnostics$pareto_k
  lwf_new <- as.vector(weights(psisf_new))

  # gather results
  list(
    lwf = lwf_new,
    kf = kf_new
  )
}




