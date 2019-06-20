library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
source('experiments/utils.R')
source('code/stanfit_funs.R')
source("code/mm_loo_funs.R")

S <- 20000
nrep <- 1
seed <- 1
batch <- 1

SEED = seed + nrep*(batch - 1)
set.seed(SEED)

raw_data = read.table('data/srrs2.dat',header = T,sep = ",")
raw_data_mn = raw_data[raw_data$state == 'PA',]
raw_data_mn$fips = raw_data_mn$stfips*1000 + raw_data_mn$cntyfips

cty = read.table('data/cty.dat',header = T,sep = ",")
cty_mn = cty[cty$st == 'PA',]
cty_mn$fips = cty_mn$stfips*1000 + cty_mn$ctfips

mn_data = merge(raw_data_mn,cty_mn[,7:8],by=c('fips'))
mn_data = mn_data[!duplicated(mn_data),]

mn_data$log_radon = log(mn_data$activity + 0.1)

n = dim(mn_data)[1]
counties = length(unique(mn_data$county))

n_county = length(unique(mn_data$county))
county = as.numeric(factor(mn_data$county,levels = unique(mn_data$county)))


standata = list(N = length(mn_data$log_radon), x = mn_data$floor, y = mn_data$log_radon,J = n_county, county = county)
stanmodel = stan_model(file='stanmodels/hier_varying_slope_and_intercept.stan')

fit <- sampling(stanmodel, data = standata, chains = 4, seed = seed, iter = S/2,cores = 1)
loo = loo(fit)














fit_rstan_mock <- stan(file = 'stanmodels/hier_varying_slope_and_intercept.stan', data = standata,
                       chains = 0)


# oracle mean and covariance
samplee = as.matrix(fit)[,1:141]
samplee[,1:3] = log(samplee[,1:3])

plot(colMeans(samplee)) # set this as mean
# similarly variance
plot(log(matrixStats::colVars(samplee)))


# gaussian proposal from advi
vbpars = c('sigma_y','sigma_a','sigma_b','a','b','mu_a','mu_b')

vb1 <- vb(stanmodel, data = standata,output_samples = 20000,pars = vbpars)
vb2 <- vb(stanmodel, data = standata,output_samples = 10000,pars = vbpars,algorithm = "fullrank")
#vb = vb1
vbsample = as.matrix(vb)[,1:141]
vbsample[,1:3] = log(vbsample[,1:3])
plot(colMeans(samplee))
points(colMeans(vbsample),col = 'red')

plot(log(matrixStats::colVars(samplee)))
points(log(matrixStats::colVars(vbsample)),col = 'red')


smean = colMeans(samplee)
cov = cov(samplee)

smean = colMeans(vbsample)
cov = cov(vbsample)

#svar = matrixStats::colVars(samplee)
#svar[svar < 0.1] = 0.1


samples <- mvtnorm::rmvnorm(4000,mean = smean,sigma = cov)


wrapper <- function(upars) {
  log(mvtnorm::dmvnorm(upars,mean = smean,sigma = cov))
}


source("rstan_code/arbitrary_proposal_alternative_stop_cond.R")
loo_rstan_arbitrary <- loo_stanfit_arbitrary(x = fit_rstan_mock, upars = samples, log_prob_prop = wrapper, k_thres = 0.7, cov = TRUE, max_iters = 10)

loo_rstan_arbitrary


#plot(loo)

loo$diagnostics$pareto_k[1:10]
loo$pointwise[1:10,1]
loo_rstan_arbitrary







# sanity test


wrapper <- function(upars) {
  log_prob_upars_stanfit(fit_rstan, upars)
}

upars =  unconstrain_pars_stanfit(fit,as.matrix(fit))


source("rstan_code/arbitrary_proposal.R")
loo_rstan_arbitrary <- loo_stanfit_arbitrary(x = fit_rstan_mock, upars =upars, log_prob_prop = wrapper, k_thres = 1, cov = TRUE)


loo_rstan_arbitrary








elpd_by_S = elpd_by_S_glm(k_threshold = 0.5,nrep = nrep,stanmodel = stanmodel,standata = standata,S_sizes = c(S),seed = SEED,control = list(adapt_delta = 0.999,max_treedepth=20))
saveRDS(object = elpd_by_S, file = paste("out/",S,"_",batch,".rds",sep=""))









test = readRDS(file='/home/topi/temp_2019_05_04/datas/radon_gauss/4000_2.rds')

test2 = readRDS(file='/home/topi/temp_2019_05_04/datas/radon/4000_2.rds')

plot(test2$psis_ks[,1,][test2$psis_ks[,1,] > 0.7], ylim = c(0,1.1))
points(test$ks[test2$psis_ks[,1,] > 0.7],col = 'red')
points(test2$ada_ks[,1,][test2$psis_ks[,1,] > 0.7],col = 'blue')


plot(test$ks)
points(test$ks_orig)



