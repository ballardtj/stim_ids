
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


##### ANALYSES ON PARAMETERS FROM PERSON LEVEL MODEL #####
library(reshape2)

load(file="data/derived/fit_vtfB_pl.RData")

parms  = rstan::extract(fit)

#1 = anodal pre, #2 = anodal post, #3 = cathodal pre, #4 = cathocal post, #5 = sham pre, #5 = sham post

v_true = melt(parms$v_true)
names(v_true) <- c('iteration','subject','condition','v_true')

B = melt(parms$B)
names(B) <- c('iteration','subject','condition','B')


hannah_dat = as.tibble(left_join(v_true,B)) %>%
  mutate(condition = factor(condition,levels=1:6,
                            labels=c('anodal.pre','anodal.dpost',
                                     'cathodal.pre','cathodal.dpost',
                                     'sham.pre','sham.dpost'))) %>%
  separate(col=condition,into=c('session','time')) %>%
  group_by(subject,session,time) %>%
  summarise_all(funs(mean=mean)) %>%
  select(-iteration_mean)

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

hannah_dat2 = left_join(hannah_dat,ids)

write_csv(hannah_dat2,path="data/clean/drift_and_thres_x_condition.csv")


regress_parms = as.tibble(left_join(v_true,B)) %>%
  mutate(condition = factor(condition,levels=1:6,
                            labels=c('anodal.pre','anodal.dpost',
                                     'cathodal.pre','cathodal.dpost',
                                     'sham.pre','sham.dpost'))) %>%
  gather(key=parm,value=value,v_true:B) %>%
  separate(col=condition,into=c('session','time')) %>%
  spread(key=time,value=value) %>%
  mutate( training_effect = pre - dpost) %>%
  select(-pre,-dpost) %>%
  spread(key=session, training_effect) %>%
  mutate(anodal_rt_disruption = sham - anodal,
         cathodal_rt_disruption = sham - cathodal) %>%
  group_by(subject,parm) %>%
  summarise(anodal_rt_disruption_m = mean(anodal_rt_disruption),
            cathodal_rt_disruption_m = mean(cathodal_rt_disruption))

regress_data_trimmed = regress_data %>%
  mutate(subject = as.numeric(as.factor(subject))) %>%
  select(subject,gender:Interaction)

regress_parms2 = left_join(regress_parms,regress_data_trimmed)


### THRESHOLD ###
#standardised, polyomial regression - cathodal
B_parms = filter(regress_parms2,parm=="B")

res=lm(cathodal_rt_disruption_m ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
         PFCwhitematter_z + PFCGlutamate_z + PFCGABA_z + PFCGlutamate_2 + PFCGABA_2 +  Interaction,data=B_parms)

summary(res)

summary(res)$coefficients[,1]

expand.grid(GABA = seq(-2,2,by=0.1),
            Glut = seq(-2,2,by=0.1)) %>%
  mutate( threshold_pred = summary(res)$coefficients[1,1] +
            summary(res)$coefficients[8,1]*Glut +
            summary(res)$coefficients[9,1]*GABA +
            summary(res)$coefficients[12,1]*Glut*GABA +
            summary(res)$coefficients[10,1]*Glut^2 +
            summary(res)$coefficients[11,1]*GABA^2) %>%
  ggplot() +
  geom_raster(aes(x=Glut,y=GABA,fill=threshold_pred)) +
  geom_contour(aes(x=Glut,y=GABA,z=threshold_pred)) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill="Threshold Disruption")


#standardised, polyomial regression - anodal
res=lm(anodal_rt_disruption_m ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
         PFCwhitematter_z + PFCGlutamate_z + PFCGABA_z + PFCGlutamate_2 + PFCGABA_2 +  Interaction,data=B_parms)

summary(res)

expand.grid(GABA = seq(-2,2,by=0.1),
            Glut = seq(-2,2,by=0.1)) %>%
  mutate( threshold_pred = summary(res)$coefficients[1,1] +
            summary(res)$coefficients[8,1]*Glut +
            summary(res)$coefficients[9,1]*GABA +
            summary(res)$coefficients[12,1]*Glut*GABA +
            summary(res)$coefficients[10,1]*Glut^2 +
            summary(res)$coefficients[11,1]*GABA^2) %>%
  ggplot() +
  geom_raster(aes(x=Glut,y=GABA,fill=threshold_pred)) +
  geom_contour(aes(x=Glut,y=GABA,z=threshold_pred)) +
  scale_fill_distiller(palette = "Spectral")  +
  labs(fill="Threshold Disruption")





### RATE ###
#standardised, polyomial regression - cathodal
v_parms = filter(regress_parms2,parm=="v_true")

res=lm(cathodal_rt_disruption_m ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
         PFCwhitematter_z + PFCGlutamate_z + PFCGABA_z + PFCGlutamate_2 + PFCGABA_2 +  Interaction,data=v_parms)

summary(res)

summary(res)$coefficients[,1]

expand.grid(GABA = seq(-2,2,by=0.1),
            Glut = seq(-2,2,by=0.1)) %>%
  mutate( rate_pred = summary(res)$coefficients[1,1] +
            summary(res)$coefficients[8,1]*Glut +
            summary(res)$coefficients[9,1]*GABA +
            summary(res)$coefficients[12,1]*Glut*GABA +
            summary(res)$coefficients[10,1]*Glut^2 +
            summary(res)$coefficients[11,1]*GABA^2) %>%
  ggplot() +
  geom_raster(aes(x=Glut,y=GABA,fill=rate_pred)) +
  geom_contour(aes(x=Glut,y=GABA,z=rate_pred)) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill="Rate Disruption")


#standardised, polyomial regression - anodal
res=lm(anodal_rt_disruption_m ~ gender_z + age_z + taskpairing_2 + taskpairing_3 + PFCgreymatter_z +
         PFCwhitematter_z + PFCGlutamate_z + PFCGABA_z + PFCGlutamate_2 + PFCGABA_2 +  Interaction,data=v_parms)

summary(res)

expand.grid(GABA = seq(-2,2,by=0.1),
            Glut = seq(-2,2,by=0.1)) %>%
  mutate( rate_pred = summary(res)$coefficients[1,1] +
            summary(res)$coefficients[8,1]*Glut +
            summary(res)$coefficients[9,1]*GABA +
            summary(res)$coefficients[12,1]*Glut*GABA +
            summary(res)$coefficients[10,1]*Glut^2 +
            summary(res)$coefficients[11,1]*GABA^2) %>%
  ggplot() +
  geom_raster(aes(x=Glut,y=GABA,fill=rate_pred)) +
  geom_contour(aes(x=Glut,y=GABA,z=rate_pred)) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill="Rate Disruption")


#Correlation between drift rate and threshold disruption

regress_parms2 %>%
  select(subject,parm,anodal_rt_disruption_m) %>%
  spread(key=parm,value=anodal_rt_disruption_m) %>%
  ggplot() +
  geom_point(aes(x=B,y=v_true))

regress_parms2 %>%
  select(subject,parm,cathodal_rt_disruption_m) %>%
  spread(key=parm,value=cathodal_rt_disruption_m) %>%
  ggplot() +
  geom_point(aes(x=B,y=v_true))


