
#clear workspace
rm(list=ls())

#load packages
library(tidyverse)
library(tidybayes)

#get files. Each one represents a participant x session combination (3 sessions per part.)
files = list.files(path="data/raw/",pattern=".txt")

#variable names
var_names=c('X1',
            'X2',
            'mapping', #1 or 2, corresponds to mapping structure. If 1, 2AFC is key 4-5. If 2, 2AFC is key 1 and 8.
            'X4',
            'X5',
            'phase', #1 is practice (before stimulation), #2 before stimulation, #3 immediate post, #4 delayed post (2 and 4 are of interest)
            'mode',#6 = 6AFC, #2 = 2AFC
            'X8',
            'trial', #trial counter for each block (all identical)
            'stimulus', #number 1 to 6 (for 2AFC, stimulus A was 1-3, stimulus B was 4-6)
            'X11',
            'response_raw', #response key (numbers 1-8, different keys for 6AFC and 2AFC)
            'correct', #1 = correct, #0 = incorrect
            'rt',  #in seconds
            'X15',
            'X16',
            'X17',
            'X18',
            'X19',
            'X20')

#For mapping 1
#stim 1 = resp 1
#stim 2 = resp 2
#stim 3 = resp 3
#stim 4 = resp 6
#stim 5 = resp 7
#stim 6 = resp 8

#For mapping 2
#stim 1 = resp 2
#stim 2 = resp 3
#stim 3 = resp 4
#stim 4 = resp 5
#stim 5 = resp 6
#stim 6 = resp 7

data_list=list()
for(i in 1:length(files)){

  tmp_dat=read.table(paste0("data/raw/",files[i]),header=T,col.names=var_names) %>%
    select(phase,mode,trial,mapping,stimulus, response_raw,correct,rt)
  strng = strsplit(files[i],split="_")[[1]]
  tmp_dat$subject = strng[3]
  tmp_dat$session = strsplit(strng[4],split=".txt")[[1]]
  data_list[[i]]=tmp_dat
}
dat = bind_rows(data_list)

imported_data_tmp = as.tibble(dat %>%
  mutate(subject = as.numeric(subject),
         session = factor(as.numeric(session),levels=1:3,labels=c('anodal','cathodal','sham')),
         phase = as.factor(phase),
         time = factor((phase=="2") +  (phase=="4")*2,levels=1:2,labels=c('pre','dpost') ),
         mode = as.factor(mode),
         mapping = as.factor(mapping),
         stimulus = as.factor(stimulus),
         response = as.factor((mapping==1)*(response_raw<=3)*(response_raw) +
                    (mapping==1)*(response_raw>3)*(response_raw-2) +
                    (mapping==2)*(response_raw-1))))


individual_diffs = read.csv("data/raw/Neurochemical_variables.csv") %>%
  mutate(subject = Subj) %>% select(-Subj)

imported_data = left_join(imported_data_tmp,individual_diffs,by="subject")

# #All subjects with 2160 trials, except 1 subject who didn't complete
# as.data.frame(
# dat %>% group_by(subject) %>%
#   summarise(count = length(trial))
# )

save(imported_data,file="data/clean/imported_data.RData")

#Trim to cut excluded subjects
#subject exluded due to poor accuracy or issues with scans (unless otherwise states)
trimmed_data_tmp = imported_data %>%
  filter(subject!=4,
         subject!=5,
         subject!=12,
         subject!=27,  #didn't finish study
         subject!=45,
         subject!=48,
         subject!=49,
         subject!=50,
         subject!=56,
         phase!="1", #cut practice trial
         phase!="3", #cut immediate post-test
         mode=="6"  #take only 6AFC
         )

trimmed_data = trimmed_data_tmp %>%
  filter(rt>0)  #negative rts indicate non-response, so cut)

#1-dim(trimmed_data)[1] / dim(trimmed_data_tmp)[1] #< 1% cut due to non-response

save(trimmed_data,file="data/clean/trimmed_data.RData")


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

COVS = cbind(
  rep(1,nrow(individual_diffs)),
  (individual_diffs$Gender==1) - (individual_diffs$Gender==2),
  scale(individual_diffs$Age),  #standardised age
  as.numeric(individual_diffs$taskpairing==2), #dummy variable representing if task pairing has value of 2
  as.numeric(individual_diffs$taskpairing==3), #dummy variable representing if task pairing has value of 3
  scale(individual_diffs$PFCgreymatter),
  scale(individual_diffs$PFCGABA),
  scale(individual_diffs$PFCGlutamate),
  scale(individual_diffs$PFCRatio),
  scale(individual_diffs$PFCGABA)*scale(individual_diffs$PFCGlutamate)
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


