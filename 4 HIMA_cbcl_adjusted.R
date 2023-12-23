#--- HIMA_hot_deck
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

# Add early-childhood CBCL total problems as a predictor
GO1momscr <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/GO1momscr.sav") 
GO2momscr <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/GO2momscr.sav") 
GO3momscr <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/GO3momscr.sav") 
GO4momscr <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/GO4momscr.sav") 

#~ Remove "le" columns (merging problems when all 4 data sets have the same set of variables, and I don't need these)
GO1momscr = GO1momscr[,!grepl("^le",names(GO1momscr))]
GO2momscr = GO2momscr[,!grepl("^le",names(GO2momscr))]
GO3momscr = GO3momscr[,!grepl("^le",names(GO3momscr))]
GO4momscr = GO4momscr[,!grepl("^le",names(GO4momscr))]

#~ Construct CES-D totals at each time point (if missing any items, missing the sum)

#~~ Merge the GO data sets
data_list <- list(GO1momscr, GO2momscr, GO3momscr, GO4momscr)
data_list_subid <- lapply(data_list, function(data) {
  data %>%
    rename(subid=famid)
})

GO <- Reduce(function(x, y) merge(x, y, by="subid", all=TRUE), 
             data_list_subid)

namesGO <- data.frame(names(GO))

#-- Create CBCL total problems for each of the 4 time points across early childhood
cbcl_early <- GO %>% 
  mutate(
    # Recode missing/refused codes in individual cbcl questions to NA
    m1cl113b = case_when(m1cl113b==999 ~ NA_real_,
                         T ~ m1cl113b)
  )

cbcl_early <- cbcl_early %>%
  mutate(
    # Sum up values (0,1,2) on individual cbcl questions at each time point
    cbcl_tot_prob1 = rowSums(.[c(257:374,376)], na.rm=T),
    cbcl_tot_prob2 = rowSums(.[c(986:1047,1049:1105,1107)], na.rm=T),
    cbcl_tot_prob3 = rowSums(.[c(1644:1708,1710:1766,1768)], na.rm=T),
    cbcl_tot_prob4 = rowSums(.[c(2397:2461,2463:2519,2521)], na.rm=T)
    ) %>%
  dplyr::select(subid,m1cl113b,cbcl_tot_prob1,cbcl_tot_prob2,cbcl_tot_prob3,cbcl_tot_prob4)

summary(cbcl_early$cbcl_tot_prob1)
summary(cbcl_early$cbcl_tot_prob2)
summary(cbcl_early$cbcl_tot_prob3)
summary(cbcl_early$cbcl_tot_prob4)

d <- merge(x=d, y=cbcl_early, all.x=T) %>%
  # Impute cbcl total problems values with median values
  mutate(
    cbcl_tot_prob1 = case_when(is.na(cbcl_tot_prob1) ~ 33,
                              T ~ cbcl_tot_prob1),
    cbcl_tot_prob2 = case_when(is.na(cbcl_tot_prob2) ~ 31,
                               T ~ cbcl_tot_prob2),
    cbcl_tot_prob3 = case_when(is.na(cbcl_tot_prob3) ~ 30,
                               T ~ cbcl_tot_prob3),
    cbcl_tot_prob4 = case_when(is.na(cbcl_tot_prob4) ~ 27,
                               T ~ cbcl_tot_prob4),
    cbcl_tot_prob = pmax(cbcl_tot_prob1, cbcl_tot_prob2, cbcl_tot_prob3, cbcl_tot_prob4)
  )

summary(d$cbcl_tot_prob)
paste0("Mean=",round(mean(d$cbcl_tot_prob),2), ", SD=", round(sqrt(var(d$cbcl_tot_prob)),2))

#----------------------------------------------------------#
#--- Standardize exposures & mediators to mean=0, std=1 ---#
#----------------------------------------------------------#

names(d)

vars<-names(d[,c(13,14,25:43,48:49)])

d <- d %>% 
  #Standardize mediators to have mean 0, std 1
  mutate_at(c(vars), ~(scale(.) %>% as.vector))


#-- Stratify data by biological sex
db <- d %>% filter(SEX==0)
dg <- d %>% filter(SEX==1)

summary(as.factor(db$TANNER_STAGE))
summary(as.factor(dg$TANNER_STAGE))


#--------------------------------#
#--- HIMA - mutually adjusted ---#
#--------------------------------#

thr_hima <- list()
for (i in names(d)[c(48,49)]){
  X = d$FIL_THREAT
  #Standardize outcome to have mean 0, std 1
  Y = d[[i]]
  print(summary(Y))
  M = d %>% dplyr::select(attentothreat,
                          ADAPTATION,STROOP_FEAR, STROOP_HAPPY, 
                          ACC_ATOM, ACC_CTOM, 
                          scr,TANNER_STAGE,
                          wasitscv, wasitscm,
                          INHIBITION_INHIBIT_BASELINE_RT, INHIBITION_SWITCH_BASELINE_RT,
                          STROOP_ACC, rs_rt, TotalStars)
  COV.XM = d %>% dplyr::select(SEX, S1AGE, POV_CHRONICITY, cesd_mom_max, cbcl_tot_prob, FIL_DEPRIVATION)
  COV.MY = d %>% dplyr::select(SEX, S1AGE, POV_CHRONICITY, FIL_DEPRIVATION, cesd_mom_max, cbcl_tot_prob)
  
  test = hima(X=X,Y=Y,M=M,COV.XM=COV.XM, COV.MY=COV.MY,
              verbose = TRUE)
  test$mediator = row.names(test)
  
  thr_hima[[i]] <- test
}

dep_hima <- list()
for (i in names(d)[c(48,49)]){
  X = d$FIL_DEPRIVATION
  #Standardize outcome to have mean 0, std 1
  Y = d[[i]]
  M = d %>% dplyr::select(attentothreat,
                          ADAPTATION,STROOP_FEAR, STROOP_HAPPY, 
                          ACC_ATOM, ACC_CTOM, 
                          scr,TANNER_STAGE,
                          wasitscv, wasitscm,
                          INHIBITION_INHIBIT_BASELINE_RT, INHIBITION_SWITCH_BASELINE_RT,
                          STROOP_ACC,rs_rt, TotalStars)
  COV.XM = d %>% dplyr::select(SEX, S1AGE, POV_CHRONICITY, cesd_mom_max, cbcl_tot_prob, FIL_THREAT)
  COV.MY = d %>% dplyr::select(SEX, S1AGE, POV_CHRONICITY, FIL_THREAT, cesd_mom_max, cbcl_tot_prob)
  
  test = hima(X=X,Y=Y,M=M,COV.XM=COV.XM, COV.MY=COV.MY,
              verbose = TRUE)
  test$mediator = row.names(test)
  
  dep_hima[[i]] <- test
}


thr<-bind_rows(thr_hima, .id = "outcome") %>%
  mutate(exposure="Threat")
dep<-bind_rows(dep_hima, .id = "outcome") %>%
  mutate(exposure = "Deprivation")
results_mutually_adj <- rbind(thr,dep) 

names(results_mutually_adj) <- c("outcome","alpha","P_alpha","beta","P_beta","gamma","gamma_SE",
                        "alpha_beta","perc_total_effect","Bonferroni.p","BH.FDR","mediator","exposure")

results_mutually_adj <- results_mutually_adj %>%
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
write.csv(results_mutually_adj, file=here("results","HIMA_results_mutually_adj_INT_EXT_TESTING_NO_NOGO_cbcl_tot_prob.csv"))

#-------------------------------------------------------#
#---            Dissect the HIMA results             ---#

noInt_est_CI<-function(coefs,row,label,type){
    coefs<-tibble::rownames_to_column(coefs, "Characteristic")
    names(coefs) <- c("Characteristic","Estimate","StdErr","T","pvalue")
    coefs<- coefs %>%
      mutate(
        sig=case_when(
          pvalue < 0.001 ~ "***",
          pvalue < 0.01 ~ "**",
          pvalue < 0.05 ~ "*",
          pvalue>=0.05 ~ "")
      )
  est<-round(coefs[row,]$Estimate,2)
  ll<-round(est-1.96*coefs[row,]$StdErr,2)
  ul<-round(est+1.96*coefs[row,]$StdErr,2)
  sig<-coefs[row,]$sig
  est_ci <- paste0(est,"(",ll,",",ul,")",sig)
  label<-label
  type<-type
  return(cbind(label=label,type=type,est_ci=est_ci))
}

#-------------------------------------------------------#
#-- OVERALL models for INT & THREAT informed by HIMA ---#

#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_thr_rs=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","gamma")
gamma_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TANNER_STAGE & rs_rt), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(rs_rt~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),4,"Reward Sensitivity (RT contrast)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),5, "Tanner Stage","beta")

THR_INT<-data.frame(rbind(gamma_thr_rs,gamma_thr_tanner, alpha_thr_rs_rt, alpha_thr_tanner, beta_thr_rs_rt, beta_thr_tanner))
THR_INT$outcome = "Latent Internalizing"
THR_INT$exposure = "Threat"
THR_INT$pop = "All"

THR_INT<-THR_INT %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 



#-------------------------------------------------------#
#-- BOYS models for INT & THREAT informed by HIMA ---#

#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_thr_rs=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","gamma")
gamma_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TANNER_STAGE & rs_rt), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(rs_rt~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),4,"Reward Sensitivity (RT contrast)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),5, "Tanner Stage","beta")

THR_INT_BOYS<-data.frame(rbind(gamma_thr_rs,gamma_thr_tanner, alpha_thr_rs_rt, alpha_thr_tanner, beta_thr_rs_rt, beta_thr_tanner))
THR_INT_BOYS$outcome = "Latent Internalizing"
THR_INT_BOYS$exposure = "Threat"
THR_INT_BOYS$pop = "Boys"

THR_INT_BOYS<-THR_INT_BOYS %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 



#-------------------------------------------------------#
#-- GIRLS models for INT & THREAT informed by HIMA ---#


#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_thr_rs=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","gamma")
gamma_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TANNER_STAGE & rs_rt), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(rs_rt~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),4,"Reward Sensitivity (RT contrast)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),5, "Tanner Stage","beta")

THR_INT_GIRLS<-data.frame(rbind(gamma_thr_rs,gamma_thr_tanner, alpha_thr_rs_rt, alpha_thr_tanner, beta_thr_rs_rt, beta_thr_tanner))
THR_INT_GIRLS$outcome = "Latent Internalizing"
THR_INT_GIRLS$exposure = "Threat"
THR_INT_GIRLS$pop = "Girls"

THR_INT_GIRLS<-THR_INT_GIRLS %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 


#------------------------------------------------------------#
#-- OVERALL models for INT & DEPRIVATION informed by HIMA ---#

#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","gamma")
gamma_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TANNER_STAGE & rs_rt), adjusted for age, sex, early-life poverty and maternal depression
alpha_rs_rt=noInt_est_CI(data.frame(summary(lm(rs_rt~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","alpha")
alpha_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),4,"Reward Sensitivity (RT contrast)", "beta")
beta_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),5, "Tanner Stage", "beta")

DEP_INT<-data.frame(rbind(gamma_rs_rt, gamma_tanner, alpha_rs_rt, alpha_tanner, beta_rs_rt, beta_tanner))
DEP_INT$outcome = "Latent Internalizing"
DEP_INT$exposure = "Deprivation"
DEP_INT$pop = "All"

DEP_INT<-DEP_INT %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 

#---------------------------------------------------------#
#-- BOYS models for INT & DEPRIVATION informed by HIMA ---#

#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","gamma")
gamma_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TANNER_STAGE & rs_rt), adjusted for age, sex, early-life poverty and maternal depression
alpha_rs_rt=noInt_est_CI(data.frame(summary(lm(rs_rt~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","alpha")
alpha_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),4,"Reward Sensitivity (RT contrast)", "beta")
beta_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),5, "Tanner Stage", "beta")

DEP_INT_BOYS<-data.frame(rbind(gamma_rs_rt, gamma_tanner, alpha_rs_rt, alpha_tanner, beta_rs_rt, beta_tanner))
DEP_INT_BOYS$outcome = "Latent Internalizing"
DEP_INT_BOYS$exposure = "Deprivation"
DEP_INT_BOYS$pop = "Boys"

DEP_INT_BOYS<-DEP_INT_BOYS %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 


#----------------------------------------------------------#
#-- GIRLS models for INT & DEPRIVATION informed by HIMA ---#

#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","gamma")
gamma_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TANNER_STAGE & rs_rt), adjusted for age, sex, early-life poverty and maternal depression
alpha_rs_rt=noInt_est_CI(data.frame(summary(lm(rs_rt~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","alpha")
alpha_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),4,"Reward Sensitivity (RT contrast)", "beta")
beta_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),5, "Tanner Stage", "beta")

DEP_INT_GIRLS<-data.frame(rbind(gamma_rs_rt, gamma_tanner, alpha_rs_rt, alpha_tanner, beta_rs_rt, beta_tanner))
DEP_INT_GIRLS$outcome = "Latent Internalizing"
DEP_INT_GIRLS$exposure = "Deprivation"
DEP_INT_GIRLS$pop = "Girls"

DEP_INT_GIRLS<-DEP_INT_GIRLS %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 

#-------------------------------------------------------#
#-- OVERALL models for EXT & THREAT informed by HIMA ---#

#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (Total stars)", "gamma")
gamma_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage", "gamma")

#-- ALPHAS: Effects of threat on retained mediators (TotalStars,TANNER_STAGE & BothRuns_All_NoGo_Trials_Accurac), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_stars=noInt_est_CI(data.frame(summary(lm(TotalStars~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),4,"Reward Sensitivity (Total stars)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),5, "Tanner Stage","beta")

THR_EXT<-data.frame(rbind(gamma_stars, gamma_tanner, alpha_thr_stars, alpha_thr_tanner, beta_thr_stars, beta_thr_tanner))
THR_EXT$outcome = "Latent Externalizing"
THR_EXT$exposure = "Threat"
THR_EXT$pop = "All"
THR_EXT<-THR_EXT %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 


#-------------------------------------------------------#
#-- BOYS models for EXT & THREAT informed by HIMA ---#

#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Reward Sensitivity (Total stars)", "gamma")
gamma_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Tanner Stage", "gamma")

#-- ALPHAS: Effects of threat on retained mediators (TotalStars,TANNER_STAGE & BothRuns_All_NoGo_Trials_Accurac), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_stars=noInt_est_CI(data.frame(summary(lm(TotalStars~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),4,"Reward Sensitivity (Total stars)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),5, "Tanner Stage","beta")

THR_EXT_BOYS<-data.frame(rbind(gamma_stars, gamma_tanner, alpha_thr_stars, alpha_thr_tanner, beta_thr_stars, beta_thr_tanner))
THR_EXT_BOYS$outcome = "Latent Externalizing"
THR_EXT_BOYS$exposure = "Threat"
THR_EXT_BOYS$pop = "Boys"

THR_EXT_BOYS<-THR_EXT_BOYS %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 

#-------------------------------------------------------#
#-- GIRLS models for EXT & THREAT informed by HIMA ---#

#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Reward Sensitivity (Total stars)", "gamma")
gamma_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Tanner Stage", "gamma")

#-- ALPHAS: Effects of threat on retained mediators (TotalStars,TANNER_STAGE & BothRuns_All_NoGo_Trials_Accurac), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_stars=noInt_est_CI(data.frame(summary(lm(TotalStars~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),4,"Reward Sensitivity (Total stars)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),5, "Tanner Stage","beta")

THR_EXT_GIRLS<-data.frame(rbind(gamma_stars, gamma_tanner, alpha_thr_stars, alpha_thr_tanner, beta_thr_stars, beta_thr_tanner))
THR_EXT_GIRLS$outcome = "Latent Externalizing"
THR_EXT_GIRLS$exposure = "Threat"
THR_EXT_GIRLS$pop = "Girls"

THR_EXT_GIRLS<-THR_EXT_GIRLS %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 

#------------------------------------------------------------#
#-- OVERALL models for EXT & DEPRIVATION informed by HIMA ---#

#-- GAMMA: Effect of deprivation on EXTT, adjusted for age, sex, early-life poverty and maternal depression
gamma_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","gamma")
gamma_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TotalStars,TANNER_STAGE & BothRuns_All_NoGo_Trials_Accurac), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_stars=noInt_est_CI(data.frame(summary(lm(TotalStars~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),4,"Reward Sensitivity (Total stars)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),5, "Tanner Stage","beta")

DEP_EXT<-data.frame(rbind(gamma_stars, gamma_tanner, alpha_thr_stars, alpha_thr_tanner, beta_thr_stars, beta_thr_tanner))
DEP_EXT$outcome = "Latent Externalizing"
DEP_EXT$exposure = "Deprivation"
DEP_EXT$pop = "All"

DEP_EXT<-DEP_EXT %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 


#---------------------------------------------------------#
#-- BOYS models for EXT & DEPRIVATION informed by HIMA ---#

#-- GAMMA: Effect of deprivation on EXTT, adjusted for age, sex, early-life poverty and maternal depression
gamma_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","gamma")
gamma_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TotalStars,TANNER_STAGE & BothRuns_All_NoGo_Trials_Accurac), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_stars=noInt_est_CI(data.frame(summary(lm(TotalStars~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),4,"Reward Sensitivity (Total stars)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=db))[["coefficients"]]),5, "Tanner Stage","beta")

DEP_EXT_BOYS<-data.frame(rbind(gamma_stars, gamma_tanner, alpha_thr_stars, alpha_thr_tanner, beta_thr_stars, beta_thr_tanner))
DEP_EXT_BOYS$outcome = "Latent Externalizing"
DEP_EXT_BOYS$exposure = "Deprivation"
DEP_EXT_BOYS$pop = "Boys"

DEP_EXT_BOYS<-DEP_EXT_BOYS %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 

#----------------------------------------------------------#
#-- GIRLS models for INT & DEPRIVATION informed by HIMA ---#

#-- GAMMA: Effect of deprivation on EXTT, adjusted for age, sex, early-life poverty and maternal depression
gamma_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","gamma")
gamma_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TotalStars,TANNER_STAGE & BothRuns_All_NoGo_Trials_Accurac), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_stars=noInt_est_CI(data.frame(summary(lm(TotalStars~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_DEPRIVATION+FIL_THREAT+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_stars=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),4,"Reward Sensitivity (Total stars)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=dg))[["coefficients"]]),5, "Tanner Stage","beta")

DEP_EXT_GIRLS<-data.frame(rbind(gamma_stars, gamma_tanner, alpha_thr_stars, alpha_thr_tanner, beta_thr_stars, beta_thr_tanner))
DEP_EXT_GIRLS$outcome = "Latent Externalizing"
DEP_EXT_GIRLS$exposure = "Deprivation"
DEP_EXT_GIRLS$pop = "Girls"

DEP_EXT_GIRLS<-DEP_EXT_GIRLS %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci) 


#-------------------------------------------------------------------
#--- Bind together all the HIMA-informed mediation results
#--- For Table 3 and Table A.4

HIMA_ALL_BY_SEX <- rbind(THR_INT, THR_INT_BOYS, THR_INT_GIRLS,
                         DEP_INT, DEP_INT_BOYS, DEP_INT_GIRLS, 
                         THR_EXT, THR_EXT_BOYS, THR_EXT_GIRLS, 
                         DEP_EXT, DEP_EXT_BOYS, DEP_EXT_GIRLS)

write.csv(HIMA_ALL_BY_SEX, file=here("/Users/Kat/Dropbox/A_DISSERTATION/Aims/Aim1/Aim1_DT/results/HIMA_ALL_BY_SEX_INT_EXT_mutually_adjusted.csv"))


########## Interactions between threat and deprivation for INT and EXT ##########
intINT<-summary(lm(INT~FIL_THREAT*FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]
intINT #p-value for interaction = 0.1210
intEXT<-summary(lm(EXT~FIL_THREAT*FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]
intEXT #p-value for interaction = 0.0105



#--------------------------------------------------------------------------------------------
#-- For Table A.3
#-- Mutually adjusted results - models for mediators adjusted for both threat and deprivation

#-- OVERALL models for INT & THREAT informed by HIMA ---#
#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_thr_rs=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","gamma")
gamma_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TANNER_STAGE & rs_rt), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(rs_rt~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (RT contrast)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),4,"Reward Sensitivity (RT contrast)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),5, "Tanner Stage","beta")

THR_INT<-data.frame(rbind(gamma_thr_rs,gamma_thr_tanner, alpha_thr_rs_rt, alpha_thr_tanner, beta_thr_rs_rt, beta_thr_tanner))
THR_INT$outcome = "Latent Internalizing"
THR_INT$exposure = "Threat"
THR_INT$pop = "All"

THR_INT_ma<-THR_INT %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci)

#-- OVERALL models for INT & DEPRIVATION informed by HIMA ---#
#-- GAMMA: Effect of deprivation on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_thr_rs=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),3,"Reward Sensitivity (RT contrast)","gamma")
gamma_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),3,"Tanner Stage","gamma")

#-- ALPHAS: Effects of deprivation on retained mediators (TANNER_STAGE & rs_rt), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(rs_rt~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),3,"Reward Sensitivity (RT contrast)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),3,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & rs_rt) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),4,"Reward Sensitivity (RT contrast)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(INT~FIL_DEPRIVATION+FIL_THREAT+rs_rt+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),5, "Tanner Stage","beta")

DEP_INT<-data.frame(rbind(gamma_thr_rs,gamma_thr_tanner, alpha_thr_rs_rt, alpha_thr_tanner, beta_thr_rs_rt, beta_thr_tanner))
DEP_INT$outcome = "Latent Internalizing"
DEP_INT$exposure = "Deprivation"
DEP_INT$pop = "All"

DEP_INT_ma<-DEP_INT %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci)


#-- OVERALL models for EXT & THREAT informed by HIMA ---#
#-- GAMMA: Effect of threat on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_thr_rs=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","gamma")
gamma_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","gamma")

#-- ALPHAS: Effects of threat on retained mediators (TANNER_STAGE & TotalStars), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(TotalStars~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Reward Sensitivity (Total stars)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),2,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & TotalStars) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),4,"Reward Sensitivity (Total stars)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),5, "Tanner Stage","beta")

THR_EXT<-data.frame(rbind(gamma_thr_rs,gamma_thr_tanner, alpha_thr_rs_rt, alpha_thr_tanner, beta_thr_rs_rt, beta_thr_tanner))
THR_EXT$outcome = "Latent Externalizing"
THR_EXT$exposure = "Threat"
THR_EXT$pop = "All"

THR_EXT_ma<-THR_EXT %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci)

#-- OVERALL models for EXT & DEPRIVATION informed by HIMA ---#
#-- GAMMA: Effect of deprivation on INT, adjusted for age, sex, early-life poverty and maternal depression
gamma_thr_rs=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),3,"Reward Sensitivity (Total stars)","gamma")
gamma_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),3,"Tanner Stage","gamma")

#-- ALPHAS: Effects of deprivation on retained mediators (TANNER_STAGE & TotalStars), adjusted for age, sex, early-life poverty and maternal depression
alpha_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(TotalStars~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),3,"Reward Sensitivity (Total stars)","alpha")
alpha_thr_tanner=noInt_est_CI(data.frame(summary(lm(TANNER_STAGE~FIL_THREAT+FIL_DEPRIVATION+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),3,"Tanner Stage","alpha")

#-- BETAs: Effects of retained mediators (TANNER_STAGE & TotalStars) on INT, adjusted for threat, deprivation, age, sex, early-life poverty and maternal depression
beta_thr_rs_rt=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),4,"Reward Sensitivity (Total stars)","beta")
beta_thr_tanner=noInt_est_CI(data.frame(summary(lm(EXT~FIL_DEPRIVATION+FIL_THREAT+TotalStars+TANNER_STAGE+SEX+S1AGE+POV_CHRONICITY+cesd_mom_max+cbcl_tot_prob, data=d))[["coefficients"]]),5, "Tanner Stage","beta")

DEP_EXT<-data.frame(rbind(gamma_thr_rs,gamma_thr_tanner, alpha_thr_rs_rt, alpha_thr_tanner, beta_thr_rs_rt, beta_thr_tanner))
DEP_EXT$outcome = "Latent Exnternalizing"
DEP_EXT$exposure = "Deprivation"
DEP_EXT$pop = "All"

DEP_EXT_ma<-DEP_EXT %>%
  dplyr::select(pop,outcome, exposure, label, type, est_ci) %>%
  pivot_wider(names_from = type, values_from = est_ci)


HIMA_ALL_mutually_adjusted <- rbind(THR_INT_ma,DEP_INT_ma, 
                                    THR_EXT_ma, DEP_EXT_ma)

write.csv(HIMA_ALL_mutually_adjusted, file=("/Users/Kat/Dropbox/A_DISSERTATION/Aims/Aim1/Aim1_DT/results/Table3_HIMA_ALL_mutually_adjusted.csv"))




#-- Co-authors insisted that the primary analysis be mutually adjusted
#-----------------------------------#
#--- HIMA - adjusted             ---#
#--- demoted from primary        ---#
#-----------------------------------#

thr_hima <- list()
#--- Run the algorithms through CBCL int/ext & YSR int/ext outcomes
for (i in names(d)[c(48,49)]){
  X = d$FIL_THREAT
  #Standardize outcome to have mean 0, std 1
  Y = d[[i]]
  print(summary(Y))
  M = d %>% dplyr::select(attentothreat,
                          ADAPTATION,STROOP_FEAR, STROOP_HAPPY, 
                          ACC_ATOM, ACC_CTOM, 
                          scr,TANNER_STAGE,
                          wasitscv, wasitscm,
                          INHIBITION_INHIBIT_BASELINE_RT, INHIBITION_SWITCH_BASELINE_RT,
                          STROOP_ACC,rs_rt, TotalStars)
  COV.XM = d %>% dplyr::select(SEX, S1AGE, POV_CHRONICITY, cesd_mom_max, cbcl_tot_prob)
  COV.MY = d %>% dplyr::select(SEX, S1AGE, POV_CHRONICITY, FIL_DEPRIVATION, cesd_mom_max, cbcl_tot_prob)
  
  test = hima(X=X,Y=Y,M=M,COV.XM=COV.XM, COV.MY=COV.MY,
              verbose = TRUE)
  test$mediator = row.names(test)
  
  thr_hima[[i]] <- test
}

dep_hima <- list()
for (i in names(d)[c(48,49)]){
  X = d$FIL_DEPRIVATION
  #Standardize outcome to have mean 0, std 1
  Y = d[[i]]
  M = d %>% dplyr::select(attentothreat,
                          ADAPTATION,STROOP_FEAR, STROOP_HAPPY, 
                          ACC_ATOM, ACC_CTOM, 
                          scr,TANNER_STAGE,
                          wasitscv, wasitscm,
                          INHIBITION_INHIBIT_BASELINE_RT, INHIBITION_SWITCH_BASELINE_RT,
                          STROOP_ACC, rs_rt, TotalStars)
  COV.XM = d %>% dplyr::select(SEX, S1AGE, POV_CHRONICITY, cesd_mom_max)
  COV.MY = d %>% dplyr::select(SEX, S1AGE, POV_CHRONICITY, FIL_THREAT, cesd_mom_max)
  
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
write.csv(results_adj, file=here("results","HIMA_results_adj_INT_EXT_TESTING_NO_NOGO_cbcl_tot_prob.csv"))

