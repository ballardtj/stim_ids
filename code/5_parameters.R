
#clear workspace
rm(list=ls())

#load packages
library(rstan)

#load fit object
load(file="data/derived/fit_vtf.RData")

#coefficients
#1) intercept
#2) gender (1 = female, -1 = male)
#3) age
#4) dummy variable task pairing = 2
#5) dummy variable task pairing = 3
#6) PFC grey matter
#7) PFC GABA
#8) PFC glutamate
#9) PFC ratio
#10) GABA X Glutamate interaction

samples = extract(fit)


COEFS_anodal_diff = samples$COEFS_anodal_true - samples$COEFS_anodal_false
COEFS_cathodal_diff = samples$COEFS_cathodal_true - samples$COEFS_cathodal_false
COEFS_sham_diff = samples$COEFS_sham_true - samples$COEFS_sham_false

#anodal
cbind(apply(COEFS_anodal_diff,2,quantile,0.025),
      apply(COEFS_anodal_diff,2,mean),
      apply(COEFS_anodal_diff,2,quantile,0.975))

#cathodal
cbind(apply(COEFS_cathodal_diff,2,quantile,0.025),
      apply(COEFS_cathodal_diff,2,mean),
      apply(COEFS_cathodal_diff,2,quantile,0.975))

#sham
cbind(apply(COEFS_sham_diff,2,quantile,0.025),
      apply(COEFS_sham_diff,2,mean),
      apply(COEFS_sham_diff,2,quantile,0.975))

#anodal v sham (diff)
cbind(apply(COEFS_anodal_diff - COEFS_sham_diff  ,2,quantile,0.025),
      apply(COEFS_anodal_diff - COEFS_sham_diff,2,mean),
      apply(COEFS_anodal_diff - COEFS_sham_diff,2,quantile,0.975))

#cathodal v sham (diffs)
cbind(apply(COEFS_cathodal_diff - COEFS_sham_diff  ,2,quantile,0.025),
      apply(COEFS_cathodal_diff - COEFS_sham_diff,2,mean),
      apply(COEFS_cathodal_diff - COEFS_sham_diff,2,quantile,0.975))

#anodal v sham (diff)
cbind(apply(samples$COEFS_anodal_true - samples$COEFS_sham_true  ,2,quantile,0.025),
      apply(samples$COEFS_anodal_true - samples$COEFS_sham_true,2,mean),
      apply(samples$COEFS_anodal_true - samples$COEFS_sham_true,2,quantile,0.975))


#cathodal v sham (true)
cbind(apply(samples$COEFS_cathodal_true - samples$COEFS_sham_true  ,2,quantile,0.025),
      apply(samples$COEFS_cathodal_true - samples$COEFS_sham_true,2,mean),
      apply(samples$COEFS_cathodal_true - samples$COEFS_sham_true,2,quantile,0.975))

apply(COEFS_cathodal_diff,2,mean)

apply(COEFS_sham_diff,2,mean)

#load stan-formatted datalist
load(file="data/clean/stan_list.RData")

#sample
fit=stan(file="models/lba_vtf.stan",
     data=stan_list,
     pars=c('A','B','tau','v_false_pre','v_true_pre','dvt_anodal','dvt_cathodal','dvt_sham'),
     #include=F,
     #cores=4
     iter=10
     )




