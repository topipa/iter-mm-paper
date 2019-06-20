#' Moment matching algorithm for updating a loo object when Pareto k estimates are high
#'
#' @param x a fitted model object
#' @param loo loo object that is modified
#' @param post_draws A function the takes \code{x} as the first argument and
#'   returns a matrix of posterior draws of the model parameters.
#' @param log_lik A function that takes \code{x} and \code{i} and returns a
#'   matrix (one columns per chain) or a vector (all chains stacked) of
#'   log-likeliood draws of the \code{i}th observation.
#' @param unconstrain_pars A function that takes arguments \code{x}, and
#'   \code{pars} and returns posterior draws on the unconstrained space based on
#'   the posterior draws on the constraint space passed via \code{pars}.
#' @param log_prob_upars A function that takes arguments \code{x} and
#'   \code{upars} and returns a matrix of log-posterior density values of the
#'   unconstrained posterior draws passed via \code{upars}.
#' @param log_lik_upars A function that takes arguments \code{x}, \code{upars},
#'   and \code{i} and returns a vector of log-likeliood draws of the \code{i}th
#'   observation based on the unconstrained posterior draws passed via
#'   \code{upars}.
#' @param max_iters Maximum number of moment matching iterations to do for each
#'   LOO fold.
#' @param k_thres Threshold value for Pareto k values above which the moment
#'   matching algorithm is used.
#' @param split Whether to do the split transform or not.
#' @param cov Also match the covariance matrix of the samples if necessary? 
#' @param cores number of cores to use for parallelization.
#' @param ... Further arguments passed to the custom functions documented above.
#' 
#' @return An updated \code{loo} object.
mm_loo_manual <- function(x, loo, post_draws, log_lik,
                          unconstrain_pars, log_prob_upars,
                          log_lik_upars, max_iters = 20,
                          k_thres = 0.5, split = FALSE,
                          cov = TRUE, cores = getOption("mc.cores", 1), 
                          ...) {
  # TODO: name mm_loo or mmloo? The latter would be similar to reloo
  # TODO: add input checks
  S <- dim(loo)[1]
  N <- dim(loo)[2]
  pars <- post_draws(x, ...)
  # transform the model parameters to unconstrained space
  upars <- unconstrain_pars(x, pars = pars, ...)
  # number of parameter in the **parameters** block only
  npars <- dim(upars)[2]
  # if more parameters than samples, do not do Cholesky transformation
  cov <- cov && S >= npars
  # compute log-probabilities of the original parameter values
  log_prob0 <- log_prob_upars(x, upars = upars, ...)
  
  # loop over all observations
  ks <- loo$diagnostics$pareto_k
  I <- which(ks > k_thres) 
  for (i in I) {
    message("Moment matching observation ", i)
    # initialize values for this LOO-fold
    uparsi <- upars
    ki <- ks[i]
    log_liki <- log_lik(x, i, ...)
    dim(log_liki) <- c(NROW(log_liki), NCOL(log_liki), 1)
    r_effi <- loo::relative_eff(exp(log_liki), cores = cores)
    dim(log_liki) <- NULL
    
    # compute log-weights per draw
    lwi <- -log_liki
    lwi <- lwi - max(lwi)
    psis_i <- SW(loo::psis(lwi, r_eff = r_effi, cores = cores))
    lwi <- as.vector(weights(psis_i))
    lwfi <- -matrixStats::logSumExp(rep(0, S))
    
    # initialize objects that keep track of the total transformation
    total_shift <- rep(0, npars)
    total_scaling <- rep(1, npars)
    if (cov) {
      total_mapping <- diag(npars)
    }
    
    iterind <- 1
    while (iterind <= max_iters && ki > k_thres) {
      # try several transformations one by one
      # if one does not work, do not apply it and try another one
      # to accept the transformation, Pareto k needs to improve
      # other possibilities: S_eff, max(lw)
      
      # 1. match means
      trans <- shift(
        x, uparsi, lwi, i = i, 
        log_prob0 = log_prob0, 
        log_prob = log_prob_upars, 
        log_lik_upars = log_lik_upars,
        r_eff = r_effi, cores = cores, ...
      )
      if (trans$ki < ki) {
        uparsi <- trans$upars
        lwi <- trans$lwi
        lwfi <- trans$lwfi
        ki <- trans$ki
        kfi <- trans$kfi
        log_liki <- trans$log_liki
        total_shift <- total_shift + trans$shift
        iterind <- iterind + 1
        next
      }
      
      # 2. match means and marginal variances
      trans <- shift_and_scale(
        x, uparsi, lwi, i = i, 
        log_prob0 = log_prob0, 
        log_prob_upars = log_prob_upars, 
        log_lik_upars = log_lik_upars,
        r_eff = r_effi, cores = cores, ...
      )
      if (trans$ki < ki) {
        uparsi <- trans$upars
        lwi <- trans$lwi
        lwfi <- trans$lwfi
        ki <- trans$ki
        kfi <- trans$kfi
        log_liki <- trans$log_liki
        total_shift <- total_shift + trans$shift
        total_scaling <- total_scaling * trans$scaling
        iterind <- iterind + 1
        next
      }
      
      # 3. match means and covariances
      if (cov) {
        trans <- shift_and_cov(
          x, uparsi, lwi, i = i, 
          log_prob0 = log_prob0, 
          log_prob_upars = log_prob_upars, 
          log_lik_upars = log_lik_upars,
          r_eff = r_effi, cores = cores, ...
        )
        if (trans$ki < ki) {
          uparsi <- trans$upars
          lwi <- trans$lwi
          lwfi <- trans$lwfi
          ki <- trans$ki
          kfi <- trans$kfi
          log_liki <- trans$log_liki
          total_shift <- total_shift + trans$shift
          total_mapping <- trans$mapping %*% total_mapping
          iterind <- iterind + 1
          next
        }
      }
      # none of the transformations improved khat
      # so there is no need to try further
      break
    }
    
    # transformations are now done
    # if we don't do split transform, or
    # if no transformations were successful
    # stop and collect values
    if (!split || (iterind == 1)) {
      elpd_loo_i <- matrixStats::logSumExp(log_liki + lwi)
    } else {
      # compute split transformation
      split_mm <- split_mm(
        x, upars, cov, total_shift, total_scaling, total_mapping, i,
        log_prob_upars = log_prob_upars, log_lik_upars = log_lik_upars,
        cores = cores, r_effi = r_effi
      )
      elpd_loo_i <- matrixStats::logSumExp(
        split_mm$log_liki_half + split_mm$lwi_half
      )
      log_liki <- split_mm$log_liki_half
      lwi <- split_mm$lwi_half
    } 
    
    
    # pointwise estimates
    # p_loo: use original p_loo, add original elpd, subtract new elpd
    loo$pointwise[i, 3] <- loo$pointwise[i, 3] + 
      loo$pointwise[i, 1] - elpd_loo_i
    # elpd_loo
    loo$pointwise[i, 1] <- elpd_loo_i
    # mcse_elpd_loo
    E_epd_i <- exp(elpd_loo_i)
    var_epd_i <- sum((exp(lwi))^2 * (exp(log_liki) - E_epd_i)^2)
    z <- rnorm(1000, mean = E_epd_i, sd = sqrt(var_epd_i))
    loo$pointwise[i, 2] <- sqrt(var(log(z[z > 0])) / r_effi)
    # looic
    loo$pointwise[i, 4] <- -2 * elpd_loo_i
    
    # diagnostics
    loo$diagnostics$pareto_k[i] <- ki
    loo$diagnostics$n_eff[i] <- 1.0 / sum(exp(2 * lwi)) * r_effi
    
    # update psis object
    if (!is.null(loo$psis_object)) {
      loo$psis_object$log_weights[, i] <- lwi
    }
  }
  # combined estimates
  loo$estimates[1, 1] <- sum(loo$pointwise[, 1])
  loo$estimates[2, 1] <- sum(loo$pointwise[, 3])
  loo$estimates[3, 1] <- sum(loo$pointwise[, 4])
  loo$estimates[1, 2] <- sqrt(N * var(loo$pointwise[, 1]))
  loo$estimates[2, 2] <- sqrt(N * var(loo$pointwise[, 3]))
  loo$estimates[3, 2] <- sqrt(N * var(loo$pointwise[, 4]))
  # loo$estimates <- loo:::table_of_estimates(loo$pointwise)
  
  # update psis object
  if (!is.null(loo$psis_object)) {
    loo$psis_object$diagnostics <- loo$diagnostics
  }
  
  # Warn if some Pareto ks are still high
  # psislw_warnings(loo$diagnostics$pareto_k)
  
  # these will be deprecated at some point
  loo$elpd_loo <- loo$estimates[1, 1]
  loo$p_loo <- loo$estimates[2, 1]
  loo$looic <- loo$estimates[3, 1]
  loo$se_elpd_loo <- loo$estimates[1, 2]
  loo$se_p_loo <- loo$estimates[2, 2]
  loo$se_looic <- loo$estimates[3, 2]
  
  loo
}


update_quantities_i <- function(x, upars, i, log_prob0, 
                                log_prob_upars, log_lik_upars,
                                r_effi, cores, ...) {
  log_prob_new <- log_prob_upars(x, upars = upars, ...)
  log_liki_new <- log_lik_upars(x, upars = upars, i = i, ...)
  # compute new raw log importance weights, unnormalized
  lwi_new <- -log_liki_new + log_prob_new - log_prob0
  lwi_new <- lwi_new - max(lwi_new)
  psis_new <- SW(loo::psis(lwi_new, r_eff = r_effi, cores = cores))
  ki_new <- psis_new$diagnostics$pareto_k
  
  lwfi_new <- log_prob_new - log_prob0
  lwfi_new <- lwfi_new - max(lwfi_new)
  psisf_new <- SW(loo::psis(lwfi_new, r_eff = r_effi, cores = cores))
  kfi_new <- psisf_new$diagnostics$pareto_k
  # pareto smoothed weights
  lwi_new <- as.vector(weights(psis_new))
  lwfi_new <- as.vector(weights(psisf_new))
  # gather results
  list(
    lwi = lwi_new, 
    lwfi = lwfi_new,
    ki = ki_new, 
    kfi = kfi_new,
    log_liki = log_liki_new
  )
}


shift <- function(x, upars, lwi, ...) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  # transform posterior draws
  upars_new <- sweep(upars, 2, shift, "+")
  # gather updated quantities
  out <- update_quantities_i(x, upars_new, ...)
  out$upars <- upars_new
  out$shift <- shift
  out
}


shift_and_scale <- function(x, upars, lwi, ...) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  mii <- exp(lwi)* upars^2
  mii <- colSums(mii) - mean_weighted^2
  scaling <- sqrt(mii / matrixStats::colVars(upars))
  # transform posterior draws
  upars_new <- sweep(upars, 2, mean_original, "-")
  upars_new <- sweep(upars_new, 2, scaling, "*")
  upars_new <- sweep(upars_new, 2, mean_weighted, "+")
  # gather updated quantities
  out <- update_quantities_i(x, upars_new, ...)
  out$upars <- upars_new
  out$shift <- shift
  out$scaling <- scaling
  out
}


shift_and_cov <- function(x, upars, lwi, ...) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  covv <- cov(upars)
  wcovv <- cov.wt(upars, wt = exp(lwi))$cov
  chol1 <- chol(wcovv)
  chol2 <- chol(covv)
  mapping <- t(chol1) %*% solve(t(chol2))
  # transform posterior draws
  upars_new <- sweep(upars, 2, mean_original, "-")
  upars_new <- tcrossprod(upars_new, mapping)
  upars_new <- sweep(upars_new, 2, mean_weighted, "+")
  colnames(upars_new) <- colnames(upars)
  # gather updated quantities
  out <- update_quantities_i(x, upars_new, ...)
  out$upars <- upars_new
  out$shift <- shift
  out$mapping <- mapping
  out
}

#' Split moment matching algorithm for doing the split transformation
#' after the moment matching algorithm has ended for a particular LOO fold
#'
#' @param x stanfit object
#' @param upars Draws from the full posterior distribution in unconstrained space
#' @param cov Also match the covariance matrix of the samples? 
#' Setting it to true may increase computation time drastically.
#' @param total_shift Accumulated translation from the moment matching
#' @param total_scaling Accumulated scaling from the moment matching
#' @param total_mapping Accumulated scaling and rotation from the moment matching
#' @param i Index for which LOO fold is being transformed
#' @param ... passed to loo::psis
#' 
#' @param r_eff Effective number of parameters
split_mm <- function(x, upars, cov, total_shift, total_scaling,
                     total_mapping, i, log_prob_upars, 
                     log_lik_upars, r_effi, cores, ...) {
  S <- dim(upars)[1]
  S_half <- as.integer(0.5 * S)
  mean_original <- colMeans(upars)
  
  # accumulated affine transformation
  upars_trans <- sweep(upars, 2, mean_original, "-")
  upars_trans <- sweep(upars_trans, 2, total_scaling, "*")
  if (cov) {
    upars_trans <- tcrossprod(upars_trans, total_mapping)
  }
  upars_trans <- sweep(upars_trans, 2, total_shift + mean_original, "+")
  
  # inverse accumulated affine transformation
  upars_trans_inv <- sweep(upars, 2, mean_original, "-")
  if (cov) {
    upars_trans_inv <- tcrossprod(upars_trans_inv, solve(total_mapping))
  }
  upars_trans_inv <- sweep(upars_trans_inv, 2, total_scaling, "/")
  upars_trans_inv <- sweep(upars_trans_inv, 2, mean_original - total_shift, "+")
  
  # first half of upars_trans_half are T(theta)
  # second half are theta
  upars_trans_half <- upars
  take <- seq_len(S_half)
  upars_trans_half[take, ] <- upars_trans[take, , drop = FALSE]
  
  # first half of upars_half_inv are theta
  # second half are T^-1 (theta)
  upars_trans_half_inv <- upars
  take <- seq_len(S)[-seq_len(S_half)]
  upars_trans_half_inv[take, ] <- upars_trans_inv[take, , drop = FALSE]
  
  log_prob_half_trans <- log_prob_upars(x, upars = upars_trans_half, ...)
  log_prob_half_trans_inv <- log_prob_upars(x, upars = upars_trans_half_inv, ...)
  log_liki_half <- log_lik_upars(x, upars = upars_trans_half, i = i, ...)
  
  # compute weights
  lwi_half <- -log_liki_half + log_prob_half_trans - 
    (log(exp(log_prob_half_trans) + 
          exp(log_prob_half_trans_inv - log(prod(total_scaling)) - 
                log(det(total_mapping)))))
  lwfi_half <- lwi_half + log_liki_half
  
  lwi_half <- lwi_half - max(lwi_half)
  lwfi_half <- lwfi_half - max(lwfi_half)
  psis_half <- SW(loo::psis(lwi_half, r_eff = r_effi, cores = cores))
  psisf_half <- SW(loo::psis(lwfi_half, r_eff = r_effi, cores = cores))
  
  # this fitting may be risky when the draws are split
  # could also fit to separate weight sets
  lwi_half <- as.vector(weights(psis_half))
  lwfi_half <- as.vector(weights(psisf_half))
  
  list(
    lwi_half = lwi_half, 
    lwfi_half = lwfi_half, 
    log_liki_half = log_liki_half
  )
}

SW <- function(expr) {
  base::suppressWarnings(expr)
}

.warn <- function(..., call. = FALSE) {
  warning(..., call. = call.)
}

.k_help <- function() {
  "See help('pareto-k-diagnostic') for details.\n"
}

psislw_warnings <- function(k) {
  if (any(k > 0.7)) {
    .warn(
      "Some Pareto k diagnostic values are too high. ",
      .k_help()
    )
  } else if (any(k > 0.5)) {
    .warn(
      "Some Pareto k diagnostic values are slightly high. ",
      .k_help()
    )
  }
}
