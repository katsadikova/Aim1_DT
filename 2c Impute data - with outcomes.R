#~~ 2c Impute data - with outcomes.R
#-----------------------------------
#~~ Impute data with INT and EXT outcomes all at once

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
library(HIMA)
library(corrplot)
library(reshape2)


# Load in the merged data file
load("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/dat.Rdata")

#~~ Candidate mediator variables: 
# -----------------------------
#   1.AB: Attention bias to threat (difference between RT_NEUTRAL_ACCURATE and RT_ANGRY_ACCURATE = attentothreat) 
#   2.ER: Emotion regulation (ADAPTATION, STROOP_FEAR, STROOP_HAPPY) 
#   3.ToM: Cognitive and affective theory of mind (ACC_ATOM, ACC_CTOM)
#   4.FC: Fear conditioning (difference between ACQ1_ASCRAmp_SR0 and ACQ1_BSCRAmp_SR0)
#   5.PT: Pubertal timing (TANNER_STAGE) 
#   6.LA: Language ability (wasitscv)
#   7.RA: Reasoning ability (wasitscm)  
#   8.IC: Inhibitory control (INHIBITION_INHIBIT_BASELINE_RT, INHIBITION_SWITCH_BASELINE_RT, STROOP_ACC, 
#                       BothRuns_All_Go_Trials_Accuracy, BothRuns_All_NoGo_Trials_Accuracy, 
#                       BothRuns_All_Accurate_Go_Trials_RT)
#   9.RS: Reward sensitivity (TotalStar, difference between RT_4star-RT_0star = rs_rt)

#~~ Create analytic data set
dat_pre_imp <-
  dat %>%
  mutate(
    thr_quart = quantcut(FINAL_THREAT, q=4),
    dep_quart = quantcut(FINAL_DEPRIVATION, q=4),
    attentothreat = RT_NEUTRAL_ACCURATE - RT_ANGRY_ACCURATE,
    rs_rt = RT_0star-RT_4star,
    scr = ACQ1_ASCRAmp_SR0 - ACQ1_BSCRAmp_SR0
  )
levels(dat_pre_imp$thr_quart) <- c(1,2,3,4)
levels(dat_pre_imp$dep_quart) <- c(1,2,3,4)


dat_pre_imp <- 
  dat_pre_imp %>%
  mutate(
    thr_low = case_when(thr_quart=='1' ~ 1,
                        TRUE ~ 0),
    thr_none = case_when(FINAL_THREAT == min(FINAL_THREAT) ~ 1,
                         TRUE ~ 0),
    dep_low = case_when(dep_quart=='1' ~ 1,
                        TRUE ~ 0)
  )

dat_pre_imp <-
  dat_pre_imp %>%
  mutate(
    #~~ Variables from the GO data sets (early childhood age 3-6)
    num_sis = case_when(m4fa004==0 ~ m4fa004,
                        m4fa004==1 ~ m4fa010),
    num_bro = case_when(m4fa005==0 ~ m4fa005,
                        m4fa005==1 ~ m4fa011),
    num_stepsis = case_when(m4fa008==0 ~ m4fa008,
                            m4fa008==1 ~ m4fa012),
    num_stepbro = case_when(m4fa009==0 ~ m4fa009,
                            m4fa009==1 ~ m4fa013),
    first_born = case_when(m1fa045 == 1 ~ 1,
                           m1fa045 > 1 ~ 0),
    early_divorce = case_when(m4fa048==0 ~ 0,
                              m4fa049==0 ~ 0,
                              m4fa049>0 ~ 1),
    preg_health = case_when(m1fa099<2 ~ 1,
                            m1fa099>1 ~ 0),
    preg_stress = case_when(m1fa100>1 ~ 1,
                            m1fa100<3 ~ 0),
    preg_alc = case_when(m1fa102==0 & m1fa103==0 ~ 0,
                         m1fa102>0 | m1fa103>0 ~1),
    preterm_lowbw = case_when(m1fa120==1 | m1fa126 == 1 ~ 1,
                              (m1fa120==0 | m1fa120==2) & m1fa126 == 0 ~ 0)
  ) %>%
  rowwise() %>%
  mutate(
    num_sib = sum(num_sis, num_bro, num_stepsis, num_stepbro, na.rm=T), 
    num_moves = sum(m1fa096, m2fa095, m3fa095, m4fa095, na.rm=T),
    preschool = sum(m1fa097,m2fa097,m3fa097,m4fa097, na.rm=T)/4,
    cesd1 = sum(m1ce001, m1ce002, m1ce003, m1ce004,
                m1ce005, m1ce006, m1ce007, m1ce008,
                m1ce009, m1ce010, m1ce011, m1ce012,
                m1ce013, m1ce014, m1ce015, m1ce016,
                m1ce017, m1ce018, m1ce019, m1ce020, na.rm = T),
    cesd2 = sum(m2ce001, m2ce002, m2ce003, m2ce004,
                m2ce005, m2ce006, m2ce007, m2ce008,
                m2ce009, m2ce010, m2ce011, m2ce012,
                m2ce013, m2ce014, m2ce015, m2ce016,
                m2ce017, m2ce018, m2ce019, m2ce020, na.rm = T),
    cesd3 = sum(m3ce001, m3ce002, m3ce003, m3ce004,
                m3ce005, m3ce006, m3ce007, m3ce008,
                m3ce009, m3ce010, m3ce011, m3ce012,
                m3ce013, m3ce014, m3ce015, m3ce016,
                m3ce017, m3ce018, m3ce019, m3ce020, na.rm = T),
    cesd4 = sum(m4ce001, m4ce002, m4ce003, m4ce004,
                m4ce005, m4ce006, m4ce007, m4ce008,
                m4ce009, m4ce010, m4ce011, m4ce012,
                m4ce013, m4ce014, m4ce015, m4ce016,
                m4ce017, m4ce018, m4ce019, m4ce020, na.rm = T)
  ) %>%
  mutate(
    cesd_mom_max = pmax(cesd1, cesd2, cesd3, cesd4),
    max_problems = pmax(YTOTAL_PROBLEMS, PTOTAL_PROBLEMS)
  )

dat_pre_imp <-
  dat_pre_imp %>%
  dplyr::select(
    subid,
    
    # Early childhood vars
    num_sib, first_born, early_divorce, num_moves, preschool,
    preg_health, preg_stress, preg_alc, preterm_lowbw, cesd_mom_max,
    
    #Exposure vars
    FINAL_THREAT, FINAL_DEPRIVATION, thr_quart, thr_none, dep_quart,
    
    # Predictors from age 11
    S1AGE, SEX, POV_CHRONICITY, INC_NEEDS, DEGREEP1P2,
    max_problems, BIO_DAD,
    
    # Mediators
    attentothreat,
    ADAPTATION,STROOP_FEAR, STROOP_HAPPY, 
    ACC_ATOM, ACC_CTOM, 
    scr,
    TANNER_STAGE,
    wasitscv, wasitscm,
    INHIBITION_INHIBIT_BASELINE_RT, INHIBITION_SWITCH_BASELINE_RT,
    STROOP_ACC,
    BothRuns_All_Go_Trials_Accuracy, BothRuns_All_NoGo_Trials_Accuracy, 
    BothRuns_All_Accurate_Go_Trials_RT, BothRuns_All_Inaccurate_NoGo_Trials_RT,
    rs_rt, TotalStars,
    
    # Outcomes
    cbcl_int_t2, cbcl_ext_t2, ysr_int_t2, ysr_ext_t2, INT, EXT)

save(dat_pre_imp, file="/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/dat_pre_imp.Rdata")
write.csv(dat_pre_imp, file="/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/dat_pre_imp.csv")
write.csv(dat_pre_imp, file="/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/DT/data/dat_pre_imp.csv")
names(dat_pre_imp)

cor(dat_pre_imp$rs_rt, dat_pre_imp$TotalStars, use="complete.obs")

#~~ Subset the mediators into a new data set to summarize
#--------------------------------------------------------
meds <- 
  dat_pre_imp %>% dplyr::select(attentothreat,
                                ADAPTATION, STROOP_FEAR, STROOP_HAPPY,
                                ACC_ATOM, ACC_CTOM,
                                scr,
                                TANNER_STAGE,
                                wasitscv,
                                wasitscm,
                                INHIBITION_INHIBIT_BASELINE_RT, INHIBITION_SWITCH_BASELINE_RT,
                                STROOP_ACC,
                                BothRuns_All_Go_Trials_Accuracy, BothRuns_All_NoGo_Trials_Accuracy, 
                                BothRuns_All_Accurate_Go_Trials_RT, BothRuns_All_Inaccurate_NoGo_Trials_RT,
                                rs_rt, TotalStars)
names(dat_pre_imp)

#~~ Look at distributions of FINAL_THREAT & FINAL_DEPRIVATION
#-----------------------------------------------------------------------------------
library(ggplot2)

t <- ggplot(dat_pre_imp, aes(x=FINAL_THREAT)) + 
  geom_density()
t
# Add mean line
t + geom_vline(aes(xintercept=mean(FINAL_THREAT)),
              color="blue", linetype="dashed", size=1)



d <- ggplot(dat_pre_imp, aes(x=FINAL_DEPRIVATION)) + 
  geom_density()
d
# Add mean line
d + geom_vline(aes(xintercept=mean(FINAL_DEPRIVATION)),
              color="blue", linetype="dashed", size=1)


#############################################################################################################################################################
#~~ IMPUTE USING MICE


#~~ Convert haven-labeled data set to a data frame that can be imputed
names(dat_pre_imp)
dat_pre_imp <- as.data.frame(dat_pre_imp)
dat_pre_imp <- sapply(dat_pre_imp, haven::zap_labels)
write.csv(dat_pre_imp, file="/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/dat_pre_imp.csv")
names(dat_pre_imp)

#~~ Characterize missingness
dim(dat_pre_imp[complete.cases(dat_pre_imp), ])
md.pattern(dat_pre_imp, plot = TRUE)
sapply(dat_pre_imp, function(x) sum(is.na(x)))

#~~ Impute everything (with outcomes)
data_imp_all <- mice(dat_pre_imp,m=20,maxit=50,meth='pmm',seed=500)
save(data_imp_all, file="/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/data_imp_all.Rdata")


###########################################################################
## ALTERNATIVELY, IMPUTE ALL VARS EXCEPT OUTCOMES, 
## AND IMPUTE OUTCOMES ONCE IN SEPARATE DATA SETS FOR EACH OUTCOME VARIABLE
#~~ Drop the outcome variables from the data to be imputed
# names(as.data.frame(dat_pre_imp))
# dat_pre_imp_nooutcomes <- dat_pre_imp[,-c(52:59)]
# dat_pre_imp_outcomes <- dat_pre_imp[,c(1,52:59)]
# summary(dat_pre_imp_outcomes)
# 
# #~~ Data imputed below and saved - if quoted out, do not rerun the imputation every time
# data_imp <- mice(dat_pre_imp_nooutcomes,m=20,maxit=50,meth='pmm',seed=500)
# save(data_imp, file="C:/Users/Kat/Documents/PHS/Research/Deprivation and Threat/Data/data_imp.Rdata")




