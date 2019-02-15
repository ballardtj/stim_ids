
#clear workspace
rm(list=ls())

#load packages
library(tidyverse)

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


regress_data = trimmed_data %>%
  group_by(subject,session,time) %>%
  summarise(rt = mean(rt)) %>%
  spread(key=time,value=rt) %>%
  mutate(rt_change = dpost-pre) %>%
  select(-pre,-dpost) %>%
  spread(key=session,rt_change) %>%
  mutate(anodal_rt_disruption = anodal - sham,
         cathodal_rt_disruption = cathodal - sham) %>%
  ungroup() %>%
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
        Interaction =  scale(ids$PFCGABA)*scale(ids$PFCGlutamate))



res=lm(anodal_rt_disruption ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
     PFCwhitematter_z + PFCGABA_z + PFCGlutamate_z + PFCGABA_2 + PFCGlutamate_2 + Interaction,data=regress_data)

summary(res)


res=lm(cathodal_rt_disruption ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
         PFCwhitematter_z + PFCGABA_z + PFCGlutamate_z + PFCGABA_2 + PFCGlutamate_2 + Interaction,data=regress_data)

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


