
#clear workspace
rm(list=ls())

#load packages
library(tidyverse)


#RUN ANALYSIS ON HANNAH'S JASP DATA
jasp_data = read_csv(file="data/clean/DataJasp.csv") %>%
  mutate(taskpairing_2 = as.numeric(taskpairing==2),
         taskpairing_3 = as.numeric(taskpairing==3))

#original specification
res=lm(CathodeVsSham ~ age + Gender + taskpairing + PFCgreymatter + PCFwhitematter +
         PFCei,data=jasp_data)

summary(res)

#task pairing as two dichotomous variables
res=lm(CathodeVsSham ~ age + Gender + taskpairing_2 + taskpairing_3 + PFCgreymatter + PCFwhitematter +
         PFCei,data=jasp_data)

summary(res)

#### REPLICATION OF THAT ANALYSIS WITH MY COMPILED DATA ####
load(file="data/clean/trimmed_data.RData")

ids = trimmed_data %>%
  group_by(subject) %>%
  summarise(gender = Gender[1],
            age = Age[1],
            taskpairing = taskpairing[1],
            PFCgreymatter = PFCgreymatter[1],
            PFCwhitematter = PFCwhitematter[1],
            PFCGABA = PFCGABA[1],
            PFCGlutamate = PFCGlutamate[1],
            PFCRatio = PFCRatio[1])

#rt_change calculated as pre - post
#discruption calculateed as sham - active
regress_data = trimmed_data %>%
  filter(correct==1) %>%
  group_by(subject,session,time) %>%
  summarise(rt = mean(rt)) %>%
  spread(key=time,value=rt) %>%
  mutate(rt_change = pre - dpost) %>%
  select(-pre,-dpost) %>%
  spread(key=session,rt_change) %>%
  mutate(anodal_rt_disruption = sham - anodal,
         cathodal_rt_disruption = sham - cathodal) %>%
  ungroup() %>%
  mutate(gender = ids$gender,
         age = ids$age,
         taskpairing_2 = as.numeric(ids$taskpairing==2),
         taskpairing_3 = as.numeric(ids$taskpairing==3),
         PFCgreymatter = ids$PFCgreymatter,
         PFCwhitematter = ids$PFCwhitematter,
         PFCei = ids$PFCGABA / ids$PFCGlutamate) %>%
  mutate(gender_z = (ids$gender==1) - (ids$gender==2),
         age_z =  scale(ids$age),  #standardised age
         taskpairing_2 =  as.numeric(ids$taskpairing==2), #dummy variable representing if task pairing has value of 2
         taskpairing_3 = as.numeric(ids$taskpairing==3), #dummy variable representing if task pairing has value of 3
         PFCgreymatter_z = scale(ids$PFCgreymatter),
         PFCwhitematter_z =  scale(ids$PFCwhitematter),
         PFCGABA_z =   scale(ids$PFCGABA),
         PFCGlutamate_z =  scale(ids$PFCGlutamate),
        PFCGABA_2 =  scale(ids$PFCGABA)^2,
        PFCGlutamate_2 =  scale(ids$PFCGlutamate)^2,
        Interaction =  scale(ids$PFCGABA)*scale(ids$PFCGlutamate),
        PFCei = scale(ids$PFCGABA / ids$PFCGlutamate))


#unstandardised, just ratio in model
res=lm(cathodal_rt_disruption ~ gender + age+ taskpairing_2 + taskpairing_3 + PFCgreymatter+
         PFCwhitematter+ PFCei ,data=regress_data)

#standardised, just ratio in model
res=lm(cathodal_rt_disruption ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
         PFCwhitematter_z + PFCei,data=regress_data)

#standardised, just gaba and glutamate
res=lm(cathodal_rt_disruption ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
         PFCwhitematter_z + PFCGlutamate_z + PFCGABA_z,data=regress_data)

summary(res)

#standardised, polyomial regression - cathodal
res=lm(cathodal_rt_disruption ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
         PFCwhitematter_z + PFCGlutamate_z + PFCGABA_z + PFCGlutamate_2 + PFCGABA_2 +  Interaction,data=regress_data)

summary(res)


expand.grid(GABA = seq(-2,2,by=0.1),
            Glut = seq(-2,2,by=0.1)) %>%
  mutate( rt_pred = 0.0299 + -0.0415288*Glut + -0.0212570*GABA + -0.0227068*Glut*GABA + 0.0141377*Glut^2 + 0.0138719*GABA^2) %>%
  ggplot() +
  geom_raster(aes(x=Glut,y=GABA,fill=rt_pred)) +
  geom_contour(aes(x=Glut,y=GABA,z=rt_pred)) +
  scale_fill_distiller(palette = "Spectral")

#standardised, polyomial regression - anodal
res=lm(anodal_rt_disruption ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
         PFCwhitematter_z + PFCGlutamate_z + PFCGABA_z + PFCGlutamate_2 + PFCGABA_2 +  Interaction,data=regress_data)

summary(res)


expand.grid(GABA = seq(-2,2,by=0.1),
            Glut = seq(-2,2,by=0.1)) %>%
  mutate( rt_pred = 0.581 + -0.028*Glut + -0.023807*GABA + -0.027106*Glut*GABA + 0.023*Glut^2 + 0.004*GABA^2) %>%
  ggplot() +
  geom_raster(aes(x=Glut,y=GABA,fill=rt_pred)) +
  geom_contour(aes(x=Glut,y=GABA,z=rt_pred)) +
  scale_fill_distiller(palette = "Spectral")


#TODO: Run unconstrained LBA to estimate parms for each person
#TODO: Bayesian models of regression




            )



summary(res)


         )


  group_by(subject) %>%
  summarise(gender = Gender[1],
            age = Age[1],
            taskpairing = taskpairing[1],
            PFCgreymatter = PFCgreymatter[1],
            PFCwhitematter = PFCwhitematter[1],
            PFCGABA = PFCGABA[1],
            PFCGlutamate = PFCGlutamate[1],
            PFCRatio = PFCRatio[1])

COVS = cbind(
  rep(1,length(unique(trimmed_data$subject))),
  (ids$gender==1) - (ids$gender==2),
  scale(ids$age),  #standardised age
  as.numeric(ids$taskpairing==2), #dummy variable representing if task pairing has value of 2
  as.numeric(ids$taskpairing==3), #dummy variable representing if task pairing has value of 3
  scale(ids$PFCgreymatter),
  scale(ids$PFCwhitematter),
  scale(ids$PFCGABA),
  scale(ids$PFCGlutamate),
  scale(ids$PFCGABA)^2,
  scale(ids$PFCGlutamate)^2,
  scale(ids$PFCGABA)*scale(ids$PFCGlutamate)
)







#### Prepare datalist for stan

#IVs
#session (3) = anodal, cathodal, sham
#time (2) = pre, dpost

#DVs - change in cathodal, change in anodal, change in sham


# vars = cbind(
#   rep(1,nrow(trimmed_data)),
#   (trimmed_data$Gender==1) - (trimmed_data$Gender==2), #switches value of 2 to -1
#   scale(trimmed_data$Age),  #standardised age
#   as.numeric(trimmed_data$taskpairing==2), #dummy variable representing if task pairing has value of 2
#   as.numeric(trimmed_data$taskpairing==3), #dummy variable representing if task pairing has value of 3
#   scale(trimmed_data$PFCgreymatter),
#   scale(trimmed_data$PFCGABA),
#   scale(trimmed_data$PFCGlutamate),
#   scale(trimmed_data$PFCRatio),
#   scale(trimmed_data$PFCGABA)*scale(trimmed_data$PFCGlutamate)
# )

ids = trimmed_data %>%
  group_by(subject) %>%
  summarise(gender = Gender[1],
            age = Age[1],
            taskpairing = taskpairing[1],
            PFCgreymatter = PFCgreymatter[1],
            PFCwhitematter = PFCwhitematter[1],
            PFCGABA = PFCGABA[1],
            PFCGlutamate = PFCGlutamate[1],
            PFCRatio = PFCRatio[1])

COVS = cbind(
  rep(1,length(unique(trimmed_data$subject))),
  (ids$gender==1) - (ids$gender==2),
  scale(ids$age),  #standardised age
  as.numeric(ids$taskpairing==2), #dummy variable representing if task pairing has value of 2
  as.numeric(ids$taskpairing==3), #dummy variable representing if task pairing has value of 3
  scale(ids$PFCgreymatter),
  scale(ids$PFCwhitematter),
  scale(ids$PFCGABA),
  scale(ids$PFCGlutamate),
  scale(ids$PFCGABA)^2,
  scale(ids$PFCGlutamate)^2,
  scale(ids$PFCGABA)*scale(ids$PFCGlutamate)
)

stan_list = list(
  Ntotal = nrow(trimmed_data),
  Nsubj = length(unique(trimmed_data$subject)),
  Nvars = ncol(COVS),
  subject = as.numeric(as.factor(trimmed_data$subject)),
  COVS = COVS,
  dpost_anodal = (trimmed_data$time=="dpost") * (trimmed_data$session=="anodal"),# * vars,
  dpost_cathodal = (trimmed_data$time=="dpost") * (trimmed_data$session=="cathodal"),# * vars,
  dpost_sham = (trimmed_data$time=="dpost") * (trimmed_data$session=="sham"),# *vars,
  RT = cbind(trimmed_data$rt,trimmed_data$correct + 1) #1 = incorrect, 2 = correct
)

save(stan_list,file="data/clean/stan_list.RData")


