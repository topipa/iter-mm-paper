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

roaches <- read.csv(file="data/roachdata.csv", header=TRUE, sep=",")
roaches$roach1 = 0.01*roaches$roach1
y = roaches$y
x = roaches[,3:5]
offset = log(roaches[,6])
n = dim(x)[1]
k = dim(x)[2]


beta_prior_scale = 2.5
alpha_prior_scale = 5.0


stanmodel = stan_model(file='stanmodels/poisson_reg.stan')
standata = list(N = n, K = k, x = as.matrix(x), y = y,offset=offset,beta_prior_scale = beta_prior_scale,alpha_prior_scale = alpha_prior_scale, loo = 0, x_loo = numeric(0), y_loo = numeric(0), offset_loo = numeric(0))


loos_exact_bf = bfloo_glm_rep(nrep = nrep,stanmodel_loo = stanmodel,standata = standata,S = S,cores=1,seed = SEED,control = list(adapt_delta = 0.95,max_treedepth=20), offset = TRUE)
saveRDS(object = loos_exact_bf, file = paste("out_naive/",S,"_",batch,".rds",sep=""))
