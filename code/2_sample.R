
#clear workspace
rm(list=ls())

#load packages
library(rstan)

#load stan-formatted datalist
load(file="data/clean/stan_list.RData")

#GABA only
stan_list$COVS = stan_list$COVS[,1:7]
stan_list$Nvars = stan_list$Nvars-3

#sample
fit=stan(file="models/lba_vtfB_reciprocal.stan",
     data=stan_list,
     pars=c('v_true','v_false','B'),
     include=F,
     cores=4,
     init_r = 1,
     #control=list(init_r=1)#adapt_delta=0.99,max_treedepth=20)
     iter=10,
     chains=1
     )

save(fit,file="data/derived/test_fit_GABA.RData")
#Check whether model works with 10 iter

#Run full analysis (2000 iter, 4 chains, 4 cores)

#See where sampling issues are occuring


