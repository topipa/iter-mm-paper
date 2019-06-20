library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
source('experiments/utils.R')
source('code/stanfit_funs.R')
source("code/mm_loo_funs.R")

S <- 2000
rep <- 1
seed <- 1
batch <- 1

SEED = 12
set.seed(SEED)

# generate toy data
data_var = 1 # std
n = as.integer(30)
x = rnorm(n=n, mean=0,sd=data_var)

xx = 1:10*2

trues = matrix(0,length(xx),rep)
bfs = matrix(0,length(xx),rep)
raws = matrix(0,length(xx),rep)
psiss = matrix(0,length(xx),rep)
adas = matrix(0,length(xx),rep)
ada50s = matrix(0,length(xx),rep)

ks_psis = matrix(0,length(xx),rep)
ks_ada = matrix(0,length(xx),rep)
kfs_ada = matrix(0,length(xx),rep)
kfs_bf = matrix(0,length(xx),rep)
ks_ada50 = matrix(0,length(xx),rep)
kfs_ada50 = matrix(0,length(xx),rep)


for (i in 1:length(xx)) {
  x[1] = xx[i]

  # analytical psoterior predictive distribution
  xmean_loo = mean(x[-1])
  xsd_loo = sd(x[-1])
  predmean_loo = xmean_loo
  predscale_loo = sqrt(1 + 1/(n-1))*xsd_loo
  preddf_loo = n-2
  lpd_loo = log(dtnew(x[1],df = preddf_loo, mean = predmean_loo, scale = predscale_loo))

  for (ii in 1:rep) {

    # full fit
    standata = list(N = n, x = x, loo = 0, x_loo = numeric(0))
    fit <- stan(file = 'stanmodels/normaldata.stan', data = standata,
                chains = 4, seed = SEED+ii, iter = as.integer(S/2),control = list(adapt_delta = 0.95))
    samples = as.data.frame(fit) # appends the chains together
    Llik = loo::extract_log_lik(fit,parameter_name="log_lik")
    looobj = loo(Llik, save_psis=TRUE)

    # raws
    lw = sumlogs_by_col(-Llik)

    # mm
    post = as.matrix(samples[,1:2])
    ada <- mm_loo_manual_ext(x = fit, k_thres = 0.5,pareto_smoothing = TRUE,cores = 1,split=FALSE)

    # smm
    ada50 <- mm_loo_manual_ext(x = fit, k_thres = 0.5,pareto_smoothing = TRUE,cores = 1,split=TRUE)

    ##############################
    # brute force

    loos_exact_bf = numeric(n)
    for (del_index in 1:1){
      standata_loo = standata
      standata_loo$N = standata$N - 1
      standata_loo$x = standata$x[-del_index]
      standata_loo$loo = 1
      standata_loo$x_loo = array(x[del_index])

      loofit <- stan(file = 'stanmodels/normaldata.stan', data = standata_loo,
                     chains = 4, seed = SEED+ii, iter = as.integer(S/2),control = list(adapt_delta = 0.95))
      loosamples = as.data.frame(loofit)
      xdens = loosamples$`loo_log_lik`
      loos_exact_bf[del_index] = matrixStats::logSumExp(xdens) - log(S)
    }


    trues[i,ii] = lpd_loo
    bfs[i,ii] = loos_exact_bf[1]
    raws[i,ii] = matrixStats::logSumExp(Llik[,1] + lw[,1])
    psiss[i,ii] = looobj$pointwise[1,1]
    adas[i,ii] = ada$pointwise[1,1]
    ada50s[i,ii] = ada50$pointwise[1,1]

    ks_psis[i,ii] = looobj$diagnostics$pareto_k[1]
    ks_ada[i,ii] = ada$diagnostics$pareto_k[1]
    kfs_ada[i,ii] = ada$diagnostics$pareto_kf[1]
    ks_ada50[i,ii] = ada50$diagnostics$pareto_k[1]
    kfs_ada50[i,ii] = ada50$diagnostics$pareto_kf[1]

    bfLlik = loo::extract_log_lik(loofit,parameter_name="loo_log_lik")
    looobj_loo = loo::loo(-bfLlik, save_psis=TRUE)
    kfs_bf[i,ii] = looobj_loo$diagnostics$pareto_k[1]


  }
}


normaldata_df = list("trues" = trues)
normaldata_df$bfs = bfs
normaldata_df$raws = raws
normaldata_df$psiss = psiss
normaldata_df$adas = adas
normaldata_df$ada50s = ada50s

normaldata_df$ks_psis = ks_psis
normaldata_df$ks_ada = ks_ada
normaldata_df$kfs_ada = kfs_ada
normaldata_df$kfs_bf = kfs_bf
normaldata_df$ks_ada50 = ks_ada50
normaldata_df$kfs_ada50 = kfs_ada50

normaldata_df$S = S

saveRDS(object = normaldata_df, file = paste("out/",S,"_",batch,".rds",sep=""))



