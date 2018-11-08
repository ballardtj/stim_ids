
#clear workspace
rm(list=ls())

#load packages
library(rstan)

#load stan-formatted datalist
load(file="data/clean/stan_list.RData")

#sample
fit=stan(file="models/lba_vtfB.stan",
     data=stan_list,
     pars=c('v_true','v_false','B'),
     include=F,
     #cores=4
     iter=10,
     chains=1
     )




