library(rstan)
library("R.matlab")
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
source('experiments/utils.R')
source('rstan_code/stanfit_funs.R')
source("rstan_code/mm_loo_funs.R")

S <- 2000
nrep <- 1
seed <- 1
batch <- 1

SEED = seed + nrep*(batch - 1)
set.seed(SEED)
ovarian <- readMat('data/ovarian.mat')
x = ovarian$x
y = ovarian$y
n = dim(x)[1]
k = dim(x)[2]


guessnumrelevcov = 20.0
slab_scale = 2.5
scale_icept = 5.0
nu_global = 1
nu_local = 1
slab_df = 1
scale_global = guessnumrelevcov/((k - guessnumrelevcov)*sqrt(n))

control = list(adapt_delta = 0.999,max_treedepth=15)

stanmodel = stan_model(file='stanmodels/horseshoe_binclas.stan')

standata = list(N = n, d = k, x = as.matrix(x), y = c(y),scale_icept = scale_icept,scale_global = scale_global,nu_global = nu_global,nu_local = nu_local,slab_scale = slab_scale, slab_df = slab_df)

elpd_by_S = elpd_by_S_glm(k_threshold = 0.7,nrep = nrep,stanmodel = stanmodel,standata = standata,S_sizes = c(S),seed = SEED,control = control)
saveRDS(object = elpd_by_S, file = paste("out/",S,"_",batch,".rds",sep=""))






