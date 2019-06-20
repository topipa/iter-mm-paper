library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
source('experiments/utils.R')
source('code/stanfit_funs.R')
source("code/mm_loo_funs.R")

S <- 2000
nrep <- 1
seed <- 1
batch <- 1

SEED = seed + nrep*(batch - 1)
set.seed(SEED)

sldata = stackloss
n = dim(sldata)[1]
k = 3
x = sldata[,1:3]
x = normalize_matrix(x)

y = sldata[,4]
y = normalize_vector(y)

beta_prior_scale = 2.5
alpha_prior_scale = 5.0

stanmodel = stan_model(file='stanmodels/linear_reg.stan')
standata = list(N = n, K = k, x = as.matrix(x), y = y,offset=offset,beta_prior_scale = beta_prior_scale,alpha_prior_scale = alpha_prior_scale, loo = 0, x_loo = numeric(0), y_loo = numeric(0))



elpd_by_S = elpd_by_S_glm(stanmodel = stanmodel,standata = standata,S_sizes = c(S),nrep = nrep,seed = SEED,k_threshold = 0.5,control = list(adapt_delta = 0.95,max_treedepth=20))
saveRDS(object = elpd_by_S, file = paste("out/",S,"_",batch,".rds",sep=""))

