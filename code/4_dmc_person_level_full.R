rm(list=ls())
library(tidyverse)

setwd("~/stim_ids")

source ("dmc/dmc.R")
load_model ("LBA","lba_B.R")

#load data
load(file="data/clean/trimmed_data.RData")

#Construct trimmed data frame which exludes missing observations
trimmed_data$S = NA
trimmed_data$S[trimmed_data$stimulus=="1"]<-"aa"
trimmed_data$S[trimmed_data$stimulus=="2"]<-"bb"
trimmed_data$S[trimmed_data$stimulus=="3"]<-"cc"
trimmed_data$S[trimmed_data$stimulus=="4"]<-"dd"
trimmed_data$S[trimmed_data$stimulus=="5"]<-"ee"
trimmed_data$S[trimmed_data$stimulus=="6"]<-"ff"
trimmed_data$S = as.factor(trimmed_data$S)
trimmed_data$R = NA
trimmed_data$R[trimmed_data$response=="1"]<-"AA"
trimmed_data$R[trimmed_data$response=="2"]<-"BB"
trimmed_data$R[trimmed_data$response=="3"]<-"CC"
trimmed_data$R[trimmed_data$response=="4"]<-"DD"
trimmed_data$R[trimmed_data$response=="5"]<-"EE"
trimmed_data$R[trimmed_data$response=="6"]<-"FF"
trimmed_data$R = as.factor(trimmed_data$R)

dmc_data = trimmed_data %>%
  filter(!is.na(R),  #remove NA values which indicate wrong key responses (27)
         phase!="3") %>% #remove second phase (phase=3) - only compare pre vs delayed post
  mutate(time = factor((phase=="2") +  (phase=="4")*2,levels=1:2,labels=c('pre','dpost') ),
         RT=rt,
         s=as.factor(as.numeric(as.factor(subject)))) %>% #Forces subject numbers to be between 0 and N
  select(s,S,R,time,session,RT,subject) %>%
  arrange(s,session,time,S)

#------------------


model <- model.dmc(p.map=list(A="1",B=c("time","session"),mean_v=c("time","session","M"),sd_v=c("time","session","M"),t0=c("time","session"), st0="1"),
   match.map=list(M=list(aa="AA",bb="BB",cc="CC",dd="DD",ee="EE",ff="FF")),
   factors=list(S=c("aa","bb","cc","dd","ee","ff"),time=c('pre','dpost'),session=c('anodal','cathodal','sham')),
   constants=c(st0=0, sd_v.dpost.sham.false = 1),
   responses=c(aa="AA",bb="BB",cc="CC",dd="DD",ee="EE",ff="FF"),
   type="norm")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"                           "B.pre.anodal"
# [3] "B.dpost.anodal"              "B.pre.cathodal"
# [5] "B.dpost.cathodal"            "B.pre.sham"
# [7] "B.dpost.sham"                "mean_v.pre.anodal.true"
# [9] "mean_v.dpost.anodal.true"    "mean_v.pre.cathodal.true"
# [11] "mean_v.dpost.cathodal.true"  "mean_v.pre.sham.true"
# [13] "mean_v.dpost.sham.true"      "mean_v.pre.anodal.false"
# [15] "mean_v.dpost.anodal.false"   "mean_v.pre.cathodal.false"
# [17] "mean_v.dpost.cathodal.false" "mean_v.pre.sham.false"
# [19] "mean_v.dpost.sham.false"     "t0.pre.anodal"
# [21] "t0.dpost.anodal"             "t0.pre.cathodal"
# [23] "t0.dpost.cathodal"           "t0.pre.sham"
# [25] "t0.dpost.sham"
#
# Constants are (see attr(,"constants") ):
#   st0 sd_v
# 0    1
#
# Model type = norm (posdrift= TRUE )


data_model <- data.model.dmc(as.data.frame(dmc_data),model)

#--------------------------------------
#Set hierarchical priors


#Note: these priors are based off the ones used in DMC tutorial 4.6. However, in that
#tutorial, there are arbitrary differences between conditions in priors (e.g., B.r1=.6 and B.r2=0.8).
#Presumably, this is because they used the priors to simulate differences in conditions
#that could later be recovered. Here, we need the priors to be the same across conditions.
#So I've fixed all threshold, mean rate and sd rates to 1.

#subject level priors
pop.mean <- c(A=1, B.pre.anodal=1, B.dpost.anodal=1,
              B.pre.cathodal=1, B.dpost.cathodal=1,
              B.pre.sham=1, B.dpost.sham=1,
              mean_v.pre.anodal.true=1, mean_v.dpost.anodal.true=1,
              mean_v.pre.cathodal.true=1, mean_v.dpost.cathodal.true=1,
              mean_v.pre.sham.true=1, mean_v.dpost.sham.true=1,
              mean_v.pre.anodal.false=1, mean_v.dpost.anodal.false=1,
              mean_v.pre.cathodal.false=1, mean_v.dpost.cathodal.false=1,
              mean_v.pre.sham.false=1, mean_v.dpost.sham.false=1,
              sd_v.pre.anodal.true=1, sd_v.dpost.anodal.true=1,
              sd_v.pre.cathodal.true=1, sd_v.dpost.cathodal.true=1,
              sd_v.pre.sham.true=1, sd_v.dpost.sham.true=1,
              sd_v.pre.anodal.false=1, sd_v.dpost.anodal.false=1,
              sd_v.pre.cathodal.false=1, sd_v.dpost.cathodal.false=1,
              sd_v.pre.sham.false=1,
              t0.pre.anodal=.3,t0.dpost.anodal=0.3,
              t0.pre.cathodal=.3,t0.dpost.cathodal=0.3,
              t0.pre.sham=.3,t0.dpost.sham=0.3)

#scales are the same as used in tutorial 4.6 (though note that they are initially specified as very small, and then multiplied by 5 on line 82)
pop.scale <- c(A=1, B.pre.anodal=1, B.dpost.anodal=1,
              B.pre.cathodal=1, B.dpost.cathodal=1,
              B.pre.sham=1, B.dpost.sham=1,
              mean_v.pre.anodal.true=1, mean_v.dpost.anodal.true=1,
              mean_v.pre.cathodal.true=1, mean_v.dpost.cathodal.true=1,
              mean_v.pre.sham.true=1, mean_v.dpost.sham.true=1,
              mean_v.pre.anodal.false=1, mean_v.dpost.anodal.false=1,
              mean_v.pre.cathodal.false=1, mean_v.dpost.cathodal.false=1,
              mean_v.pre.sham.false=1, mean_v.dpost.sham.false=1,
              sd_v.pre.anodal.true=1, sd_v.dpost.anodal.true=1,
              sd_v.pre.cathodal.true=1, sd_v.dpost.cathodal.true=1,
              sd_v.pre.sham.true=1, sd_v.dpost.sham.true=1,
              sd_v.pre.anodal.false=1, sd_v.dpost.anodal.false=1,
              sd_v.pre.cathodal.false=1, sd_v.dpost.cathodal.false=1,
              sd_v.pre.sham.false=1,
              t0.pre.anodal=.3,t0.dpost.anodal=0.3,
              t0.pre.cathodal=.3,t0.dpost.cathodal=0.3,
              t0.pre.sham=.3,t0.dpost.sham=0.3)

p.prior <- prior.p.dmc(
  dists = rep("tnorm",36),
  p1=pop.mean,p2=pop.scale,
  lower=c(0,0,0,0,0,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,0,0,0,0,0,0,0,0,0,0,0,.1,.1,.1,.1,.1,.1),
  upper=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,1,1,1,1,1)
)

#population level priors
# mu.prior <- prior.p.dmc(
#   dists = rep("tnorm",20),
#   p1=pop.mean,
#   p2=c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,1), #mean rate distributions have an sd of 2. All others have sd of 1.
#   lower=c(0,0,0,0,0,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,.1),
#   upper=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1)
# )

# sigma.prior <- prior.p.dmc(
#   dists = rep("beta", length(p.prior)),
#   p1=c(A=1, B.pre.anodal=1, B.dpost.anodal=1,
#                 B.pre.cathodal=1, B.dpost.cathodal=1,
#                 B.pre.sham=1, B.dpost.sham=1,
#                 mean_v.pre.anodal.true=1, mean_v.dpost.anodal.true=1,
#                 mean_v.pre.cathodal.true=1, mean_v.dpost.cathodal.true=1,
#                 mean_v.pre.sham.true=1, mean_v.dpost.sham.true=1,
#                 mean_v.pre.anodal.false=1, mean_v.dpost.anodal.false=1,
#                 mean_v.pre.cathodal.false=1, mean_v.dpost.cathodal.false=1,
#                 mean_v.pre.sham.false=1, mean_v.dpost.sham.false=1,t0=1),
#   p2=rep(1,20) #All sd priors are uniform bounded between 0 and 1.
# )

# sigma.prior <- prior.p.dmc(
#   dists = rep("tnorm", length(p.prior)),
#   p1=c(A=0, B.pre.anodal=0, B.dpost.anodal=0,
#        B.pre.cathodal=0, B.dpost.cathodal=0,
#        B.pre.sham=0, B.dpost.sham=0,
#        mean_v.pre.anodal.true=0, mean_v.dpost.anodal.true=0,
#        mean_v.pre.cathodal.true=0, mean_v.dpost.cathodal.true=0,
#        mean_v.pre.sham.true=0, mean_v.dpost.sham.true=0,
#        mean_v.pre.anodal.false=0, mean_v.dpost.anodal.false=0,
#        mean_v.pre.cathodal.false=0, mean_v.dpost.cathodal.false=0,
#        mean_v.pre.sham.false=0, mean_v.dpost.sham.false=0,t0=0),
#   p2=rep(1,20), #All sd priors are unit half normals
#   lower=rep(0,20),
#   upper=rep(NA,20)
# )


# Make a hyper-prior list
# pp.prior <- list(mu.prior, sigma.prior)

#Save priors
# prior_list=list(p.prior=p.prior,pp.prior=pp.prior)

# ------------------------------------------------
# Get starting values

starting_samples <- h.samples.dmc(nmc = 100, p.prior,data_model, thin = 5)#, pp.prior = pp.prior)
#save(starting_samples, file = "2_hierarchical_starting_samples_v2.RData")

unstuck_samples <- h.run.unstuck.dmc(starting_samples, p.migrate = .05, cores = 7)
#save(unstuck_samples, file = "2_hierarchical_unstuck_samples_v2.RData")

final_samples <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=unstuck_samples), nmc=100,cores=7,finalrun=T,finalI=500,minN=400,meanN=450)
save(final_samples, file = "data/derived/dmc_final_samples_person_level_full.RData")


#----------------------------------------#
#             CONVERGENCE                #
#----------------------------------------#

gelman.diag.dmc(final_samples)

# 14   40    6   26   18   20   30   25   38   19   33   32    7   13   36   43   15
# 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.09 1.10 1.10 1.10
# 1   44   51   55   31   22   17   21   11    2   52    3   34   16   39   42   28
# 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10
# 54   23   53   35   29   37    9   46   10    8   41   47   24
# 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10
# Mean
# [1] 1.1

effectiveSize.dmc(final_samples)

#most are well above 1000

plot.dmc(final_samples,subject=47)


#-------------------------------------------------#
#             POSTERIOR PREDITIVES                #
#-------------------------------------------------#

pp=h.post.predict.dmc(samples=final_samples,save.simulation=T,cores=7)

sim = do.call(rbind, pp) %>%
  mutate(prev_reps = lag(reps),
         new_subject = as.numeric(reps < prev_reps))

sim$new_subject[1] = 1
sim$subject = cumsum(sim$new_subject)

### Accuracy ###

sim = sim %>%
  #mutate(s = rownames(sim),
 #        s = gsub("\\..*","",s)) %>%
  mutate(correct = as.numeric(S == tolower(R))) %>%
  group_by(subject,time,session,reps) %>%
  mutate(accuracy = mean(correct)) %>%
  group_by(time,session,reps) %>%
  summarise(mean_accuracy = mean(accuracy)) %>%
  group_by(time,session) %>%
  summarise(prop_m = mean(mean_accuracy),
            prop_l = quantile(mean_accuracy,0.025),
            prop_u = quantile(mean_accuracy,0.975),
            source = "Model")

data = lapply(pp, function(x) attr(x, "data"))
data = do.call(rbind, data)
data = data %>%
  #mutate(s = rownames(data),
  #       s = gsub("\\..*","",s)) %>%
  mutate(correct = as.numeric(S == tolower(R))) %>%
  group_by(subject,time,session) %>%
  mutate(accuracy = mean(correct)) %>%
  group_by(time,session) %>%
  summarise(mean_accuracy = mean(accuracy)) %>%
  group_by(time,session) %>%
  summarise(prop_m = mean(mean_accuracy),
            prop_l = NA,#prop_m - sd(prop)/sqrt(n()),
            prop_u = NA,#prop_m + sd(prop)/sqrt(n()),
            source = "Data")

pp_smry =  bind_rows(data,sim) %>%
  ungroup() #%>%
  # mutate(Coh = factor(Coh,levels=c('c25','c20','c15','c10'),labels=c('Very\nEasy','Easy','Hard','Very\nHard')),
  #        St = factor(St,levels=c('left','right'),labels=c('Leftward Motion','Rightward Motion')),
  #        Emph = factor(Emph,levels=c('speed','accuracy'),labels=c('Speed Emphasis','Accuracy Emphasis')))

pp_accuracy = ggplot(data=pp_smry,aes(x=time,y=prop_m,group=source,colour=source)) +
  geom_errorbar(aes(ymax = prop_u, ymin = prop_l), width= 0.2) +
  geom_point(pch=21, size=2) +
  geom_line(aes(group=source)) +
  ylab("Proportion Correct") + xlab('Time') +
  scale_y_continuous(breaks = seq(0.5,1,0.1),limits = c(0.5,1)) +
  facet_grid(~session)






#+ theme_minimal()





h.pick.stuck.dmc(final_samples)

#



ps=summary.dmc(final_samples)
write.csv(ps[[1]]$statistics,"hierarchical_parameters.csv")

pp=h.post.predict.dmc(samples=final_samples)

plot.pp.dmc(pp,style="cdf")


#Variable transformations - mean parameters
Nsubj=length(ps)
Nvar=length(ps[[1]]$statistics[,'Mean'])+1 #extra for subject
subject_means=data.frame(matrix(NA,Nsubj,Nvar))
names(subject_means)=c('subject',names(ps[[1]]$statistics[,'Mean']))
for(i in 1:Nsubj){
  subject_means[i,1] = as.numeric(names(ps[i]))
  subject_means[i,2:Nvar]=ps[[i]]$statistics[,'Mean']
}

#Difference in Drift
subject_means$mean_v.pre.anodal.diff = subject_means$mean_v.pre.anodal.true      - subject_means$mean_v.pre.anodal.false
subject_means$mean_v.dpost.anodal.diff = subject_means$mean_v.dpost.anodal.true  - subject_means$mean_v.dpost.anodal.false
subject_means$mean_v.pre.cathodal.diff = subject_means$mean_v.pre.cathodal.true  - subject_means$mean_v.pre.cathodal.false
subject_means$mean_v.dpost.cathodal.diff = subject_means$mean_v.dpost.cathodal.true - subject_means$mean_v.dpost.cathodal.false
subject_means$mean_v.pre.sham.diff = subject_means$mean_v.pre.sham.true - subject_means$mean_v.pre.sham.false
subject_means$mean_v.dpost.sham.diff = subject_means$mean_v.dpost.sham.true - subject_means$mean_v.dpost.sham.false

#Sum of Drifts
subject_means$mean_v.pre.anodal.sum = subject_means$mean_v.pre.anodal.true      + subject_means$mean_v.pre.anodal.false
subject_means$mean_v.dpost.anodal.sum = subject_means$mean_v.dpost.anodal.true  + subject_means$mean_v.dpost.anodal.false
subject_means$mean_v.pre.cathodal.sum = subject_means$mean_v.pre.cathodal.true  + subject_means$mean_v.pre.cathodal.false
subject_means$mean_v.dpost.cathodal.sum = subject_means$mean_v.dpost.cathodal.true + subject_means$mean_v.dpost.cathodal.false
subject_means$mean_v.pre.sham.sum = subject_means$mean_v.pre.sham.true + subject_means$mean_v.pre.sham.false
subject_means$mean_v.dpost.sham.sum = subject_means$mean_v.dpost.sham.true + subject_means$mean_v.dpost.sham.false

#Change in Drift from Pre to Delayed Post
subject_means$mean_v.change.anodal.diff = subject_means$mean_v.dpost.anodal.diff - subject_means$mean_v.pre.anodal.diff
subject_means$mean_v.change.cathodal.diff = subject_means$mean_v.dpost.cathodal.diff - subject_means$mean_v.pre.cathodal.diff
subject_means$mean_v.change.sham.diff = subject_means$mean_v.dpost.sham.diff - subject_means$mean_v.pre.sham.diff

subject_means$mean_v.change.anodal.sum = subject_means$mean_v.dpost.anodal.sum - subject_means$mean_v.pre.anodal.sum
subject_means$mean_v.change.cathodal.sum = subject_means$mean_v.dpost.cathodal.sum - subject_means$mean_v.pre.cathodal.sum
subject_means$mean_v.change.sham.sum = subject_means$mean_v.dpost.sham.sum - subject_means$mean_v.pre.sham.sum

subject_means$mean_v.change.anodal.true = subject_means$mean_v.dpost.anodal.true - subject_means$mean_v.pre.anodal.true
subject_means$mean_v.change.cathodal.true = subject_means$mean_v.dpost.cathodal.true - subject_means$mean_v.pre.cathodal.true
subject_means$mean_v.change.sham.true = subject_means$mean_v.dpost.sham.true - subject_means$mean_v.pre.sham.true

subject_means$mean_v.change.anodal.false = subject_means$mean_v.dpost.anodal.false - subject_means$mean_v.pre.anodal.false
subject_means$mean_v.change.cathodal.false = subject_means$mean_v.dpost.cathodal.false - subject_means$mean_v.pre.cathodal.false
subject_means$mean_v.change.sham.false = subject_means$mean_v.dpost.sham.false - subject_means$mean_v.pre.sham.false

subject_means$B.change.anodal = subject_means$B.dpost.anodal - subject_means$B.pre.anodal
subject_means$B.change.cathodal = subject_means$B.dpost.cathodal - subject_means$B.pre.cathodal
subject_means$B.change.sham = subject_means$B.dpost.sham - subject_means$B.pre.sham

write.csv(subject_means,"subject_means.csv",row.names=F)

subject_means %>%
  summarise_all("mean")


x = trimmed_data %>%
  group_by(subject) %>%
  summarise(rt = mean(rt)) %>%
  ungroup() %>%
  mutate(newsubject = as.factor(as.numeric(subject)))



s=as.factor(as.numeric(subject))) %>% #Forces subject numbers to be between 0 and N
  select(s,S,R,time,session,RT) %>%
  arrange(s,session,time,S)


 filter(trial==1,phase==2,session=="anodal")


#
#
# #final_samples <- h.run.dmc(samples.dmc(nmc=300, samples=burnin_samples,thin=5),cores=7)
# #save(final_samples, file = "../model_output/1_group_level_final_samples.RData")
# #load( file = "../model_output/1_group_level_final_samples.RData")
#
#
#
# plot.dmc(final_samples,pll.chain=TRUE,start=1)
#
# gelman.diag.dmc(final_samples)
#
# pick.stuck.dmc(final_samples)
#
# effectiveSize.dmc(final_samples)
#
# pp=post.predict.dmc(samples=final_samples,style)
# plot.pp.dmc(pp,style="cdf")
#
# summary.dmc(final_samples)
#
#
