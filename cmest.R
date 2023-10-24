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
library(CMAverse)
library(here)

#-- Set working directory
setwd("/Users/Kat/Dropbox/A_DISSERTATION/Aims/Aim1/Aim1_DT")

# Hot-deck imputed data
d<-data.frame(read_sas("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/DT/data/dat_hotdeck_int_ext.sas7bdat"))
summary(as.factor(d$ImpIndex)) #141/227 observations imputed due to missingness in 1 or more variables
names(d)


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
names(d)

#----------------------------------------------------------#
#--- Standardize exposures & mediators to mean=0, std=1 ---#
#----------------------------------------------------------#
vars<-names(d[,c(13,14,25:43)])

d <- d %>% 
  #Standardize mediators to have mean 0, std 1
  mutate_at(c(vars), ~(scale(.) %>% as.vector))

#-- Stratify data by biological sex
db <- d %>% filter(SEX==0)
dg <- d %>% filter(SEX==1)

#-- Checking normality of mediator and outcomes to assess appropriateness of using linear regressions
hist(d$rs_rt) # skewed - long left tail - but fine
hist(d$INT) # fine

#--- Regression estimated mediation - not adjusted for Tanner
mediation.rb <- cmest(data = d, 
                      model = "rb", 
                      outcome = names(d)[48], 
                      exposure = "FIL_THREAT", 
                      mediator = c("rs_rt"),  
                      EMint = F,
                      basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max","cbcl_tot_prob"), 
                      mreg = list("linear"), 
                      yreg = "linear", 
                      a = 1, 
                      astar = 0, 
                      mval = list(1),
                      estimation = "imputation", 
                      inference = "bootstrap")
summary(mediation.rb)$summarydf

#--- Regression estimated mediation -> additionally adjusting for Tanner - THIS IS WHAT's REPORTED IN FIG 1!!!!!
mediation.rb <- cmest(data = d, 
                      model = "rb", 
                      outcome = names(d)[48], 
                      exposure = "FIL_THREAT", 
                      mediator = c("rs_rt"),  
                      EMint = F,
                      basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max","cbcl_tot_prob","TANNER_STAGE"), 
                      mreg = list("linear"), 
                      yreg = "linear", 
                      a = 1, 
                      astar = 0, 
                      mval = list(1),
                      estimation = "imputation", 
                      inference = "bootstrap")
summary(mediation.rb)$summarydf


#--- G-formula estimated mediation with Tanner allowed to be affected by threat (tho it's not in this sample) - no substantial change
mediation.gform <- cmest(data = d, 
                         model = "gformula", 
                         outcome = names(d)[48], 
                         exposure = "FIL_THREAT", 
                         mediator = c("rs_rt"),  
                         EMint = F,
                         basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max","cbcl_tot_prob"), 
                         postc = c("TANNER_STAGE"),
                         mreg = list("linear"), 
                         yreg = "linear", 
                         postcreg = list("linear"),
                         a = 1, 
                         astar = 0, 
                         mval = list(1),
                         estimation = "imputation", 
                         inference = "bootstrap")
summary(mediation.gform)$summarydf

#-----------------------------------------------------------------------------------#
#-----------------  Testing other things/sensitivity analyses  ---------------------#

#--- Regression estimated mediation, adjusted for deprivation & Tanner
mediation.rb.d <- cmest(data = d, 
                      model = "rb", 
                      outcome = names(d)[48], 
                      exposure = "FIL_THREAT", 
                      mediator = c("rs_rt"),  
                      EMint = F,
                      basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max", "cbcl_tot_prob", "FIL_DEPRIVATION"),
                      mreg = list("linear"), 
                      yreg = "linear", 
                      a = 1, 
                      astar = 0, 
                      mval = list(1),
                      estimation = "imputation", 
                      inference = "bootstrap")
summary(mediation.rb.d)$summarydf

#--- Regression estimated mediation - check for CBCL internalizing
mediation.rb.cbcl <- cmest(data = d, 
                      model = "rb", 
                      outcome = names(d)[44], 
                      exposure = "FIL_THREAT", 
                      mediator = c("rs_rt"),  
                      EMint = F,
                      basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max","cbcl_tot_prob"), 
                      mreg = list("linear"), 
                      yreg = "linear", 
                      a = 1, 
                      astar = 0, 
                      mval = list(1),
                      estimation = "imputation", 
                      inference = "bootstrap")
summary(mediation.rb.cbcl)$summarydf

#--- Regression estimated mediation - check for YSR internalizing
mediation.rb.ysr <- cmest(data = d, 
                           model = "rb", 
                           outcome = names(d)[46], 
                           exposure = "FIL_THREAT", 
                           mediator = c("rs_rt"),  
                           EMint = F,
                           basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max","cbcl_tot_prob"), 
                           mreg = list("linear"), 
                           yreg = "linear", 
                           a = 1, 
                           astar = 0, 
                           mval = list(1),
                           estimation = "imputation", 
                           inference = "bootstrap")
summary(mediation.rb.ysr)$summarydf



#--- Among boys only - not enough power really
#--- G-formula estimated mediation with baseline max_problems allowed to be affected by threat
mediation.gform <- cmest(data = db, 
                         model = "gformula", 
                         outcome = names(d)[48], 
                         exposure = "FIL_THREAT", 
                         mediator = c("rs_rt"),  
                         EMint = F,
                         basec = c("S1AGE","POV_CHRONICITY","cesd_mom_max"), 
                         postc = c("TANNER_STAGE"),
                         mreg = list("linear"), 
                         yreg = "linear", 
                         postcreg = list("linear"),
                         a = 1, 
                         astar = 0, 
                         mval = list(1),
                         estimation = "imputation", 
                         inference = "bootstrap")
summary(mediation.gform)$summarydf

#--- Among girls only - not enough power really
#--- G-formula estimated mediation with baseline max_problems allowed to be affected by threat
mediation.gform <- cmest(data = dg, 
                         model = "gformula", 
                         outcome = names(d)[48], 
                         exposure = "FIL_THREAT", 
                         mediator = c("rs_rt"),  
                         EMint = F,
                         basec = c("S1AGE","POV_CHRONICITY","cesd_mom_max"), 
                         postc = c("TANNER_STAGE"),
                         mreg = list("linear"), 
                         yreg = "linear", 
                         postcreg = list("linear"),
                         a = 1, 
                         astar = 0, 
                         mval = list(1),
                         estimation = "imputation", 
                         inference = "bootstrap")
summary(mediation.gform)$summarydf 

summary(lm(data=test1, zEXT~lp_lowExec+lp_lowToMscr))

tab <- test1 %>% group_by(lp_cat) %>%
  summarize(mean_t=mean(zFINAL_THREAT),
            mean_d=mean(zFINAL_DEPRIVATION))

cor.test(test1$lp_lowToMscr,test1$zFINAL_THREAT)
