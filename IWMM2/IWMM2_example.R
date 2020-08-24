
source('IWMM2/IWMM2.R')

# todo check if something is unconstrained

set.seed(5)

S <- 4000

prop_mean <- c(0,0)
prop_var <- c(0.1,0.1)

target_mean <- c(5,5)
target_var <- c(1,1)


prop_sample <- mvtnorm::rmvnorm(S,mean = prop_mean, sigma = diag(prop_var))
prop_density <- function(upars, ...) {
  mvtnorm::dmvnorm(upars, mean = prop_mean, sigma = diag(prop_var), log = TRUE)
}

target_density <- function(upars ,...) {
  mvtnorm::dmvnorm(upars, mean = target_mean, sigma = diag(target_var), log = TRUE)
}



iw_mean <- IW2(prop_sample, expectation_fun =  function(upars, ...) {
  upars[,1:2]},
  log_prob_prop_upars_fun = prop_density,
  log_prob_target_upars_fun = target_density)




iw2_mean <- IW2(prop_sample, expectation_fun =  function(upars, ...) {
  upars[,1:2]},
  log_prob_prop_upars_fun = prop_density,
  log_prob_target_upars_fun = target_density,
  moment_match = TRUE, k_threshold = 0.3)




# true means are (5,5)
iw_mean
iw2_mean





iw_moment2 <- IW2(prop_sample, expectation_fun =  function(upars, ...) {
  upars[,1:2]^2},
  log_prob_prop_upars_fun = prop_density,
  log_prob_target_upars_fun = target_density)




iw2_moment2 <- IW2(prop_sample, expectation_fun =  function(upars, ...) {
  upars[,1:2]^2},
  log_prob_prop_upars_fun = prop_density,
  log_prob_target_upars_fun = target_density,
  moment_match = TRUE, k_threshold = 0.3)


iw_var <- iw_moment2$expectation - iw_mean$expectation^2
iw2_var <- iw2_moment2$expectation - iw2_mean$expectation^2

# true variances are (1,1)
iw_var
iw2_var



##





