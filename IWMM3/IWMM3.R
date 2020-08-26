#' Generic importance weighted moment matching (IWMM) importance sampling algorithm.
#' Can compute arbitrary expectations when we have a sample from the distribution
#' over which the expectation is defined
#'
#'
#' @param upars A matrix of model parameters in unconstrained space.
#' @param lw_psis A psis object from the loo package computed from importance weights
#' @param expectation_fun A function whose expectation is being computed.
#' The function takes arguments `upars`.
#' @param log_prob_upars_fun
#' @param log_expectation_fun Logical indicating whether the expectation_fun
#' returns its values as logarithms or not. Defaults to FALSE. If set to TRUE,
#' the expectation function must be nonnegative (before taking the logarithm).
#' @param k_threshold Threshold for pareto k values above which the moment matching is used.
#' Set to 0.5 by default. Lower (higher) threshold leads
#' to better (worse) accuracy but slower (faster) computation.
#' @param cov_transform Logical; Indicate whether to match the covariance matrix of the
#'   samples or not. If `FALSE`, only the mean and marginal variances are
#'   matched.
#' @param dim which dimension the expectation depends on. not used at the moment
#' @param ... Further arguments passed to the custom functions documented above.
#'
#' @return Returns a list with 3 elements: the expectation of interest,
#' and 2 pareto k diagnostic values.
IWMM3 <- function(upars, log_prob_upars_fun, expectation_fun,
                  log_expectation_fun = FALSE, k_threshold = 0.5, cov_transform = TRUE, ...) {
  orig_log_prob <- log_prob_upars_fun(upars = upars, ...)
  npars <- ncol(upars)
  S <- nrow(upars)
  cov_transform <- cov_transform && S >= 10 * npars

  lw <- rep(0,nrow(upars))
  lw_orig <- lw



  upars2 <- upars
  total_shift2 <- rep(0, npars)
  total_scaling2 <- rep(1, npars)
  total_mapping2 <- diag(npars)



  if (log_expectation_fun) {
    lwf <- lw + expectation_fun(upars2, ...)
    # Here we expect that expectation_fun is nonnegative (before log)')
  }
  else {
    lwf <- lw + log(abs(expectation_fun(upars2, ...)))
  }
  psisf <- suppressWarnings(loo::psis(lwf))
  kf <- psisf$diagnostics$pareto_k


  if (ncol(lwf) > 1) {
    stop('This function is not yet supported for expectation functions that return a matrix.
         As a workaround, you can wrap your function call using apply.')
  }

  lwf <- as.vector(weights(psisf))

  while (kf > k_threshold) {

    # 1. match means
    trans <- shift(upars2, lwf)
    quantities <- update_quantities3b(
      upars2 = trans$upars,
      orig_log_prob = orig_log_prob,
      log_prob_upars_fun = log_prob_upars_fun,
      expectation_fun = expectation_fun,
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
    trans <- shift_and_scale(upars2, lwf)
    quantities <- update_quantities3b(
      upars2 = trans$upars,
      orig_log_prob = orig_log_prob,
      log_prob_upars_fun = log_prob_upars_fun,
      expectation_fun = expectation_fun,
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

    if (cov_transform) {
      # 3. match means and covariances
      trans <- shift_and_cov(upars2, lwf)
      quantities <- update_quantities3b(
        upars2 = trans$upars,
        orig_log_prob = orig_log_prob,
        log_prob_upars_fun = log_prob_upars_fun,
        expectation_fun = expectation_fun,
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
    }



    break
  }

  S_half <- as.integer(0.5 * S)


  # original parameters
  mean_original <- colMeans(upars)

  # accumulated affine transformation
  upars_T2 <- sweep(upars, 2, mean_original, "-")
  upars_T2 <- sweep(upars_T2, 2, total_scaling2, "*")
  if (cov_transform) {
    upars_T2 <- tcrossprod(upars_T2, total_mapping2)
  }
  upars_T2 <- sweep(upars_T2, 2, total_shift2 + mean_original, "+")

  # inverse accumulated affine transformation
  upars_T2inv <- sweep(upars, 2, mean_original, "-")
  if (cov_transform) {
    upars_T2inv <- tcrossprod(upars_T2inv, solve(total_mapping2))
  }
  upars_T2inv <- sweep(upars_T2inv, 2, total_scaling2, "/")
  upars_T2inv <- sweep(upars_T2inv, 2, mean_original - total_shift2, "+")




  # these are the real used draws
  # first half of upars_trans are theta
  # second half are T2(theta)
  upars_trans <- upars_T2
  take <- seq_len(S_half)
  upars_trans[take, ] <- upars[take, , drop = FALSE]


  # then we need one set of pseudo draws
  # inverse transformed upars_trans
  # first half are T2^-1(theta)
  # second half are theta
  upars_trans_inv2 <- upars
  upars_trans_inv2[take, ] <- upars_T2inv[take, , drop = FALSE]


  log_prob_trans <- log_prob_upars_fun(upars = upars_trans, ...)
  log_prob_trans_inv2 <- log_prob_upars_fun(upars = upars_trans_inv2, ...)


  lw_trans <-  log_prob_trans -
    log(
      exp(log_prob_trans) +
      exp(log_prob_trans_inv2 - log(prod(total_scaling2))  - log(det(total_mapping2)))
    )


  lw_trans_psis <- suppressWarnings(loo::psis(lw_trans))
  lw_trans <- as.vector(weights(lw_trans_psis))


  # replace upars and lw
  lw <- lw_trans
  upars <- upars_trans



  if (log_expectation_fun) {
    expectation <- exp(matrixStats::colLogSumExps(lw + expectation_fun(upars, ...)))
  }
  else {
    w <- exp(lw)
    expectation <- colSums(w * expectation_fun(upars, ...))
  }

  list("expectation" = expectation, "pareto_kf" = kf)
}






####################### HELPER FUNCTIONS



shift <- function(upars, lw, ...) {
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

shift_and_scale <- function(upars, lw, ...) {
  # compute moments using log weights
  S <- dim(upars)[1]
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lw) * upars)
  shift <- mean_weighted - mean_original
  mii <- exp(lw)* upars^2
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

shift_and_cov <- function(upars, lw, ...) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lw) * upars)
  shift <- mean_weighted - mean_original
  covv <- stats::cov(upars)
  wcovv <- stats::cov.wt(upars, wt = exp(lw))$cov
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


update_quantities3b <- function(upars2, orig_log_prob,
                               log_prob_upars_fun,
                               expectation_fun, log_expectation_fun,
                               ...) {


  log_prob_new <- log_prob_upars_fun(upars = upars2, ...)
  lw_new <- log_prob_new - orig_log_prob

  if (log_expectation_fun) {
    lwf_new <- lw_new + expectation_fun(upars2, ...)
  }
  else {
    lwf_new <- lw_new + log(abs(expectation_fun(upars2, ...)))
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




