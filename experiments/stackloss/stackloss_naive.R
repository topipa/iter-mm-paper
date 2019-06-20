library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
source('experiments/utils.R')

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

loos_exact_bf = bfloo_glm_rep(nrep = nrep,stanmodel_loo = stanmodel,standata = standata,S = S,cores=1,seed = SEED,control = list(adapt_delta = 0.95,max_treedepth=20))
saveRDS(object = loos_exact_bf, file = paste("out_naive/",S,"_",batch,".rds",sep=""))

