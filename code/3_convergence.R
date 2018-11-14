
#clear workspace
rm(list=ls())

#load packages
library(rstan)

##### GABA Model #####

#load fit object
load(file="data/derived/fit_vtfB_gaba.RData")

fit

#uniform priors on hyper sds

# Inference for Stan model: lba_vtfB.
# 4 chains, each with iter=2000; warmup=1000; thin=1;
# post-warmup draws per chain=1000, total post-warmup draws=4000.
#
# mean se_mean    sd    2.5%     25%     50%     75%  97.5% n_eff Rhat
# A_mean                       0.56    0.00  0.16    0.22    0.47    0.58    0.67   0.85  1923 1.00
# A_sd                         0.74    0.00  0.15    0.44    0.62    0.74    0.86   0.98  2496 1.00
# B_pre_mean                   1.11    0.00  0.17    0.76    1.00    1.12    1.23   1.42  2013 1.00
# B_pre_sd                     0.44    0.00  0.17    0.19    0.32    0.41    0.52   0.86  2040 1.00
# tau_mean                     0.11    0.00  0.01    0.10    0.10    0.11    0.12   0.14  1321 1.00
# tau_sd                       0.02    0.00  0.01    0.00    0.01    0.01    0.02   0.05   517 1.00
# v_false_pre_mean             0.86    0.01  0.18    0.51    0.75    0.86    0.98   1.20   490 1.01
# v_false_pre_sd               0.42    0.00  0.15    0.19    0.31    0.40    0.50   0.79  1760 1.00
# v_true_pre_mean              3.37    0.00  0.09    3.20    3.31    3.37    3.43   3.55  3223 1.00
# v_true_pre_sd                0.11    0.00  0.09    0.00    0.04    0.09    0.16   0.35  1636 1.00

#changed to unit normals

#if fails, try centering regression coefs?


sd_pars = names(fit)[grep("_sd",names(fit))]

traceplot(fit,pars=sd_pars )

samples = extract(fit)

iter = 1

dB = samples$dB[1,]

dB_min = samples$dB_min[1]

pre_B = samples$B_pre[1,]

actual_pre_B = pre_B - dB_min

post_B = actual_pre_B + dB


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




