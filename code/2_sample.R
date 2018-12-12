
#clear workspace
rm(list=ls())

#load packages
library(rstan)

#load stan-formatted datalist
# load(file="data/clean/stan_list_few_subjects.RData")
#
# #GABA only
# stan_list$COVS = stan_list$COVS[,1:7]
# stan_list$Nvars = stan_list$Nvars-3
#
# #sample
# fit=stan(file="models/lba_vtfB_reciprocal.stan",
#      data=stan_list,
#      pars=c('v_true','v_false','B'),
#      include=F,
#      cores=7,
#      #init_r = 1,
#      control=list(max_treedepth=20),
#      #iter=10,
#      chains=7
#      )
#
# save(fit,file="data/derived/fit_vtfB_gaba_reciprocal.RData")

#Glut only
# load(file="data/clean/stan_list.RData")
# stan_list$COVS = stan_list$COVS[,c(1:6,8)]
# stan_list$Nvars = stan_list$Nvars-3
#
# #sample
# fit=stan(file="models/lba_vtfB_reciprocal.stan",
#          data=stan_list,
#          pars=c('v_true','v_false','B'),
#          include=F,
#          cores=7,
#          #init_r = 1,
#          control=list(max_treedepth=20),
#          #iter=10,
#          chains=7
# )
#
# save(fit,file="data/derived/fit_vtfB_glut_reciprocal.RData")

#Ratio only
load(file="data/clean/stan_list_few_subjects.RData")
stan_list$COVS = stan_list$COVS[,c(1:6,9)]
stan_list$Nvars = stan_list$Nvars-3

#sample
fit=stan(file="models/lba_vtfB_reciprocal_vf_fixed.stan",
         data=stan_list,
         pars=c('v_true','v_false','B'),
         include=F,
         cores=4,
         #init_r = 1,
         control=list(max_treedepth=20),
         iter=10,
         chains=4
)

save(fit,file="data/derived/fit_vtfB_ratio_reciprocal.RData")



#Check whether model works with 10 iter

#Run full analysis (2000 iter, 4 chains, 4 cores)

#See where sampling issues are occuring


