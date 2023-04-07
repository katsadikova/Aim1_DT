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

#--- Regression estimated mediation -> THIS IS WHAT's REPORTED!!!!
mediation.rb <- cmest(data = d, 
                      model = "rb", 
                      outcome = names(d)[48], 
                      exposure = "FIL_THREAT", 
                      mediator = c("rs_rt"),  
                      EMint = F,
                      basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max"), 
                      mreg = list("linear"), 
                      yreg = "linear", 
                      a = 1, 
                      astar = 0, 
                      mval = list(1),
                      estimation = "imputation", 
                      inference = "bootstrap")
summary(mediation.rb)$summarydf

#-----------------------------------------------------------------------------------#
#-----------------  Testing other things/sensitivity analyses  ---------------------#

#--- Regression estimated mediation, adjusted for deprivation & Tanner
mediation.rb.d <- cmest(data = d, 
                      model = "rb", 
                      outcome = names(d)[48], 
                      exposure = "FIL_THREAT", 
                      mediator = c("rs_rt"),  
                      EMint = F,
                      basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max", "FIL_DEPRIVATION"),
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
                      basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max"), 
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
                           basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max"), 
                           mreg = list("linear"), 
                           yreg = "linear", 
                           a = 1, 
                           astar = 0, 
                           mval = list(1),
                           estimation = "imputation", 
                           inference = "bootstrap")
summary(mediation.rb.ysr)$summarydf


#--- G-formula estimated mediation with baseline max_problems allowed to be affected by threat
mediation.gform <- cmest(data = d, 
                          model = "gformula", 
                          outcome = names(d)[48], 
                          exposure = "FIL_THREAT", 
                          mediator = c("rs_rt"),  
                          EMint = F,
                          basec = c("S1AGE","SEX","POV_CHRONICITY","cesd_mom_max"), 
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
