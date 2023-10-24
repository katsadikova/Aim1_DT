#--- HIMA_BOYS.R
#--- Date: 2/10/2023
#--------------------------------------------------
#~~ Altered HIMA function with SE output
source("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/NEW_CODE/TESTING_HIMA - output with standard errors.R")
#~~ Function to calculate SEs from estimate and p-value (https://rdrr.io/github/MathiasHarrer/dmetar/src/R/SE_from_p.R)
source("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/NEW_CODE/se.from.p.R")

library(tidyverse)
library(kableExtra)
library(gtsummary)
library(expss)
library(haven)
library(sjlabelled)
library(readxl)
library(gtools)
library(tableone)
library(mice)
#library(HIMA)
library(corrplot)
library(reshape2)
library(dplyr)
library(here)

#--- Set working directory
setwd("/Users/Kat/Dropbox/A_DISSERTATION/Aims/Aim1/Aim1_DT/")

# Hot-deck imputed data
d<-data.frame(read_sas("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/DT/data/dat_hotdeck_int_ext.sas7bdat"))
summary(as.factor(d$ImpIndex)) #141/227 observations imputed due to missingness in 1 or more variables

str(d)
names(d)

summary(d$SEX) #1=Female
#----------------------------------------------------------#
#--- Standardize exposures & mediators to mean=0, std=1 ---#
#----------------------------------------------------------#

vars<-names(d[,c(13,14,25:43,48:49)])

d <- d %>%
  #Standardize mediators to have mean 0, std 1
  mutate_at(c(vars), ~(scale(.) %>% as.vector))


#-- Stratify data by biological sex
db <- d %>% filter(SEX==0)
dg <- d %>% filter(SEX==1)

summary(as.factor(db$TANNER_STAGE))
summary(as.factor(dg$TANNER_STAGE))

#-----------------------#
#--- HIMA - adjusted ---#
#-----------------------#

thr_hima <- list()
#--- Run the algorithms through CBCL int/ext & YSR int/ext outcomes
for (i in names(db)[c(48,49)]){
  X = db$FIL_THREAT
  #Standardize outcome to have mean 0, std 1
  Y = db[[i]]
  print(summary(Y))
  M = db %>% dplyr::select(attentothreat,
                          ADAPTATION,STROOP_FEAR, STROOP_HAPPY, 
                          ACC_ATOM, ACC_CTOM, 
                          scr,TANNER_STAGE,
                          wasitscv, wasitscm,
                          INHIBITION_INHIBIT_BASELINE_RT, INHIBITION_SWITCH_BASELINE_RT,
                          STROOP_ACC,rs_rt, TotalStars)
  COV.XM = db %>% dplyr::select(S1AGE, POV_CHRONICITY, cesd_mom_max)
  COV.MY = db %>% dplyr::select(S1AGE, POV_CHRONICITY, FIL_DEPRIVATION, cesd_mom_max)
  
  test = hima(X=X,Y=Y,M=M,COV.XM=COV.XM, COV.MY=COV.MY,
              verbose = TRUE)
  test$mediator = row.names(test)
  
  thr_hima[[i]] <- test
}

dep_hima <- list()
for (i in names(db)[c(48,49)]){
  X = db$FIL_DEPRIVATION
  #Standardize outcome to have mean 0, std 1
  Y = db[[i]]
  M = db %>% dplyr::select(attentothreat,
                          ADAPTATION,STROOP_FEAR, STROOP_HAPPY, 
                          ACC_ATOM, ACC_CTOM, 
                          scr,TANNER_STAGE,
                          wasitscv, wasitscm,
                          INHIBITION_INHIBIT_BASELINE_RT, INHIBITION_SWITCH_BASELINE_RT,
                          STROOP_ACC, rs_rt, TotalStars)
  COV.XM = db %>% dplyr::select(S1AGE, POV_CHRONICITY, cesd_mom_max)
  COV.MY = db %>% dplyr::select(S1AGE, POV_CHRONICITY, FIL_THREAT, cesd_mom_max)
  
  test = hima(X=X,Y=Y,M=M,COV.XM=COV.XM, COV.MY=COV.MY,
              verbose = TRUE)
  test$mediator = row.names(test)
  
  dep_hima[[i]] <- test
}


thr<-bind_rows(thr_hima, .id = "outcome") %>%
  mutate(exposure="Threat")
dep<-bind_rows(dep_hima, .id = "outcome") %>%
  mutate(exposure = "Deprivation")
results_adj <- rbind(thr,dep) 

names(results_adj) <- c("outcome","alpha","P_alpha","beta","P_beta","gamma","gamma_SE",
                        "alpha_beta","perc_total_effect","Bonferroni.p","BH.FDR","mediator","exposure")

results_adj <- results_adj %>%
  mutate(
    alpha_LLCI = se.from.p(alpha, P_alpha, 227, effect.size.type = 'difference', calculate.g = FALSE)$LLCI,
    alpha_ULCI = se.from.p(alpha, P_alpha, 227, effect.size.type = 'difference', calculate.g = FALSE)$ULCI,
    beta_LLCI = se.from.p(beta, P_beta, 227, effect.size.type = 'difference', calculate.g = FALSE)$LLCI,
    beta_ULCI = se.from.p(beta, P_beta, 227, effect.size.type = 'difference', calculate.g = FALSE)$ULCI,
    gamma_LLCI = gamma - 1.96*gamma_SE,
    gamma_ULCI = gamma + 1.96*gamma_SE,
    alpha_CI = paste0(round(alpha,2),"(",round(alpha_LLCI,2),",", round(alpha_ULCI,2),")"),
    beta_CI = paste0(round(beta,2),"(",round(beta_LLCI,2),",", round(beta_ULCI,2),")"),
    gamma_CI = paste0(round(gamma,2),"(",round(gamma_LLCI,2),",", round(gamma_ULCI,2),")")
  ) %>%
  group_by(outcome,exposure) %>%
  mutate(
    n_meds=n()
  ) %>%
  ungroup() %>%
  mutate(
    unadj_p = Bonferroni.p/n_meds,
    joint_sig = case_when(unadj_p<0.1 ~ 1,
                          T ~ 0)) %>%
  dplyr::select(c(outcome, exposure, mediator, n_meds, alpha_CI, 
                  beta_CI, gamma_CI, alpha_beta,
                  perc_total_effect,BH.FDR,Bonferroni.p,unadj_p,joint_sig))
write.csv(results_adj, file=here("/Users/Kat/Dropbox/A_DISSERTATION/Aims/Aim1/Aim1_DT/results/HIMA_BOYS_results_adj_INT_EXT_TESTING_NO_NOGO.csv"))
