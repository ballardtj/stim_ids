rm(list=ls())
#setwd("~/Fingers/Exp1a/dmc_hierarchical/")
source ("dmc/dmc.R")
load_model ("LBA","lba_B.R")

#load data
load("data/clean/trimmed_data.RData")

#Construct trimmed data frame which exludes missing observations
dmc_data = trimmed_data %>%
 mutate(R = factor(correct,levels=c(1,0),labels=c('COR','INC')),
        s = as.factor(as.numeric(subject)), #this forces subject numbers to be between 1 and N.
        S = factor(rep(1,n()),levels=c(1,0),c('cor','inc')),
        RT = rt,
        TIME = factor(time,levels=c('pre','dpost')),
        SESSION = factor(session,levels=c('anodal','cathodal','sham'))) %>%
 select(s,S,TIME,SESSION,R,RT) %>%
 arrange(s,S)


#------------------
#Set up model

 model <- model.dmc(p.map=list(A="1",B=c("TIME","SESSION"),mean_v=c("TIME","SESSION", "M"), sd_v="1",t0=c("TIME","SESSION"), st0="1"),
                   match.map=list(M=list(cor="COR",inc="INC")),
                   factors=list(S=c("cor","inc"),TIME=c('pre','dpost'),SESSION=c('anodal','cathodal','sham')),
                   constants=c(st0=0, mean_v.pre.sham.false = 1), # ANY negative value => Inf RT/never wins
                   responses=c(correct="COR",incorrect="INC"),
                   type="norm")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"                           "B.pre.anodal"                "B.dpost.anodal"
# [4] "B.pre.cathodal"              "B.dpost.cathodal"            "B.pre.sham"
# [7] "B.dpost.sham"                "mean_v.pre.anodal.true"      "mean_v.dpost.anodal.true"
# [10] "mean_v.pre.cathodal.true"    "mean_v.dpost.cathodal.true"  "mean_v.pre.sham.true"
# [13] "mean_v.dpost.sham.true"      "mean_v.pre.anodal.false"     "mean_v.dpost.anodal.false"
# [16] "mean_v.pre.cathodal.false"   "mean_v.dpost.cathodal.false" "mean_v.dpost.sham.false"
# [19] "sd_v"                        "t0.pre.anodal"               "t0.dpost.anodal"
# [22] "t0.pre.cathodal"             "t0.dpost.cathodal"           "t0.pre.sham"
# [25] "t0.dpost.sham"
#
# Constants are (see attr(,"constants") ):
#   st0 mean_v.pre.sham.false
# 0                     1

#Model type = norm (posdrift= TRUE )

data_model <- data.model.dmc(as.data.frame(dmc_data),model)

# #--------------------------------------
# #Setpriors

#Note: these priors are based off the ones used in DMC tutorial 4.6. However, in that
#tutorial, there are arbitrary differences between conditions in priors (e.g., B.r1=.6 and B.r2=0.8).
#Presumably, this is because they used the priors to simulate differences in conditions
#that could later be recovered. Here, we need the priors to be the same across conditions.
#So I've fixed all threshold, mean rate and sd rates to 1.

#subject level priors
prior.mean <- c(A=1,
              B.pre.anodal=1, B.dpost.anodal=1, B.pre.cathodal=1, B.dpost.cathodal=1, B.pre.sham=1, B.dpost.sham=1,
              mean_v.pre.anodal.true=1, mean_v.dpost.anodal.true=1, mean_v.pre.cathodal.true=1,
              mean_v.dpost.cathodal.true=1, mean_v.pre.sham.true=1, mean_v.dpost.sham.true=1,
              mean_v.pre.anodal.false=1, mean_v.dpost.anodal.false=1, mean_v.pre.cathodal.false=1,
              mean_v.dpost.cathodal.false=1, mean_v.dpost.sham.false=1,
              sd_v = 1,
              t0.pre.anodal=0.3, t0.dpost.anodal=0.3, t0.pre.cathodal=0.3,
              t0.dpost.cathodal=0.3, t0.pre.sham=0.3, t0.dpost.sham=0.3)

pop.scale <- c(A=1,
              B.pre.anodal=1, B.dpost.anodal=1, B.pre.cathodal=1, B.dpost.cathodal=1, B.pre.sham=1, B.dpost.sham=1,
              mean_v.pre.anodal.true=1, mean_v.dpost.anodal.true=1, mean_v.pre.cathodal.true=1,
              mean_v.dpost.cathodal.true=1, mean_v.pre.sham.true=1, mean_v.dpost.sham.true=1,
              mean_v.pre.anodal.false=1, mean_v.dpost.anodal.false=1, mean_v.pre.cathodal.false=1,
              mean_v.dpost.cathodal.false=1, mean_v.dpost.sham.false=1,
              sd_v = 1,
              t0.pre.anodal=0.5, t0.dpost.anodal=0.5, t0.pre.cathodal=0.5,
              t0.dpost.cathodal=0.5, t0.pre.sham=0.5, t0.dpost.sham=0.5)


p.prior <- prior.p.dmc(
 dists = rep("tnorm",length(pop.mean)),
 p1=pop.mean,p2=pop.scale,
 lower=c(0,0,0,0,0,0,0, #A and Bs
         rep(NA,11), #mean_vs
         0, #sd_v
         rep(0.1,6)), #t0
 upper=c(rep(NA,7), #A and Bs
         rep(NA,11), #mean_vs
         NA, #sd_v
         rep(1,6)) #t0
)

#population level priors
# mu.prior <- prior.p.dmc(
#  dists = rep("tnorm",11),
#  p1=pop.mean,
#  p2=c(1,1,1,2,2,2,2,1,1,1,1), #mean rate distributions have an sd of 2. All others have sd of 1.
#  lower=c(0,0,0,NA,NA,NA,NA,0,0,0,.1),upper=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1)
# )
#
# sigma.prior <- prior.p.dmc(
#  dists = rep("beta", length(p.prior)),
#  p1=c(A=1, B.MM=1, B.NM=1,
#       mean_v.mm.true=1, mean_v.nm.true=1, mean_v.mm.false=1, mean_v.nm.false=1,
#       sd_v.mm.true = 1, sd_v.nm.true = 1, sd_v.mm.false = 1,
#       t0=1),p2=rep(1,11) #All sd priors are uniform bounded between 0 and 1.
# )
#
#
# # Make a hyper-prior list
# pp.prior <- list(mu.prior, sigma.prior)

#Save priors
# prior_list=list(p.prior=p.prior,pp.prior=pp.prior)
# save(prior_list,file="Exp1a_prior_list.RData")

# # ------------------------------------------------
# # Get starting values
#
starting_samples <- h.samples.dmc(nmc = 100, p.prior,data_model, thin = 10)

unstuck_samples <- h.run.unstuck.dmc(starting_samples, p.migrate = .05, cores = 7)

final_samples <- h.run.converge.dmc(samples=unstuck_samples, cores=7, verbose=TRUE, report=10, finalrun=TRUE,finalI=500, addtofinal = FALSE)

