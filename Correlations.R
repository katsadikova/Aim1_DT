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
library(miceadds)
library(HIMA)
library(corrplot)
library(reshape2)
library(TAM)
library(ggplot2)
library(gridExtra)
library(ggpubr) # Combine figures
library(here)
library(Hmisc)

#--- Set working directory
setwd("/Users/Kat/Dropbox/A_DISSERTATION/Aims/Aim1/Aim1_DT")


#--- Load in pre-imputation data
load("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/dat_pre_imp.Rdata")
names(dat_pre_imp)

#--- Load in hot-deck imputed data
d<-data.frame(read_sas("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/DT/data/dat_hotdeck_int_ext.sas7bdat"))
summary(as.factor(d$ImpIndex)) #141/227 observations imputed due to missingness in 1 or more variables
#-----------------------------------------------------------------------
#--- Distributions of variables and missing data counts (Appendix table)

#--- Create dummy variable for those who experienced poverty ever during early life 

dat_pre_imp <- dat_pre_imp %>%
  mutate(
    pov_any = case_when(POV_CHRONICITY>=1 ~ 1,
                        POV_CHRONICITY==0 ~ 0,
                        T ~ NA_real_),
    thr_any = 1-thr_none
  )

table(as.factor(dat_pre_imp$POV_CHRONICITY), as.factor(dat_pre_imp$pov_any))
table(as.factor(dat_pre_imp$thr_any),as.factor(dat_pre_imp$thr_none))

names(dat_pre_imp)
vars <- names(dat_pre_imp[,c(17,18,19,49,11,12,50,13,24:36,41,42,47,48)])
catVars <- names(dat_pre_imp[,c(18,49,50)])
vars
catVars
by_thr <- CreateTableOne(vars = vars, factorVars = catVars,
                         addOverall = FALSE,
                         data = dat_pre_imp)
by_thr <- data.frame(print(by_thr, missing=T))
by_thr <- tibble::rownames_to_column(by_thr, "Characteristic")

# Column labels
names(by_thr) <- c("Characteristic", "Overall", "Missing")
by_thr

# Row labels
table_rows <- c("n",  
                "Age, baseline, mean(SD)",
                "Female biological sex, n(%)",
                "Chronicity of poverty, early childhood, mean(SD)",
                "Ever poverty, early childhood, n(%)",
                "Maternal depression, early childhood, mean(SD)",
                "Threat, mean(SD)",
                "Any threat, n(%)",
                "Deprivation, mean(SD)",
                "AB: Attention bias threat, mean(SD)",
                "ER: Adaptation to emotional conflict, mean(SD)",
                "ER: Stroop - fear, mean(SD)",
                "ER: Stroop - happy, mean(SD)",
                "ToM: Accuracy on affective trials, mean(SD)",
                "ToM: Accuracy on cognitive trials, mean(SD)",
                "FC: Skin conductance response to CS+ vs CS-, mean(SD)",
                "PT: Tanner stage, mean(SD)",
                "AL: Language ability, mean(SD)",
                "AR: Reasoning ability, mean(SD)",
                "IC: Reaction time on inhibit trials, mean(SD)",
                "IC: Reaction time on switch trials, mean(SD)",
                "IC: Accuracy on Stroop task, mean(SD)",
                "RS: Reaction time on high- vs low-reward trials, mean(SD)",
                "RS: Total stars, mean(SD)",
                "Internalizing symptoms, mean(SD)",
                "Externalizing symptoms, mean(SD)")
by_thr$Characteristic <- table_rows
kbl(by_thr,
    caption = "Distributions of child characteristics, with % missing pre-imputation") %>% 
  kable_styling()

write.csv(by_thr, file=here("results","Distributions_with_missing.csv"))

#~~ Quick look at normality of outcomes
hist(d$cbcl_int_t2) #Good as is, slightly better with log transformation
hist(d$cbcl_ext_t2) #Medium - but will do
hist(d$ysr_int_t2) #Good as is
hist(d$ysr_ext_t2) #Good as is
hist(d$INT) #Good as is
hist(d$EXT) #Good as is

names(d)
#~~ Look at zero-order correlations between DT & outcomes
DT_outcomes_corr <- d[,c(13,14,48:49)]
round(cor(DT_outcomes_corr),2)[,c(1,2)]


#~~ Look at distributions of mediators & zero-order correlations between DT & mediators
DT_mediators_corr <- d[,c(13,14,25:37,42,43,48:49)]
names(DT_mediators_corr)
dim(DT_mediators_corr)

# Only between meds and DT, outcomes
med_corrs1 <- as.data.frame(round(cor(DT_mediators_corr,method = c("pearson")),2)[c(1:19),c(1,2,18,19)])

# p-values
med_corrs1_pvals <- as.data.frame(round(rcorr(as.matrix(DT_mediators_corr),type = "pearson")[["P"]],4)[c(1:19),c(1,2,18,19)])

colnames(med_corrs1) <- c("Threat", "Deprivation",
                         "Latent Internalizing", "Latent Externalizing")
rownames(med_corrs1) <- c( "Threat",
                           "Deprivation",
                           "AB: Attention bias threat",
                           "ER: Adaptation to emotional conflict",
                           "ER: Stroop - fear",
                           "ER: Stroop - happy",
                           "ToM: Accuracy on affective trials",
                           "ToM: Accuracy on cognitive trials",
                           "FC: Skin conductance response to CS+ vs CS-",
                           "PT: Tanner stage",
                           "AL: Language ability",
                           "AR: Reasoning ability",
                           "IC: Reaction time on inhibit trials",
                           "IC: Reaction time on switch trials",
                           "IC: Accuracy on Stroop task",
                           "RS: Reaction time on high- vs low-reward trials",
                           "RS: Total stars",
                           "Latent Internalizing",
                           "Latent Externalizing")

kbl(med_corrs1,
    caption = "Correlations between candidate mediators and deprivation/threat exposures and adolescent psychiatric outcomes") %>% 
  kable_styling()

write.csv(med_corrs1, file=here("results","med_corrs_INT_EXT.csv"))


# Corrs between mediators
med_corrs2 <- as.data.frame(round(cor(DT_mediators_corr),2)[c(3:10,13:15,11,12,16,17),c(3:10,13:15,11,12,16,17)])
# p-values
med_corrs2_pvals <- as.data.frame(round(rcorr(as.matrix(DT_mediators_corr),type = "pearson")[["P"]],4)[c(3:10,13:15,11,12,16,17),c(3:10,13:15,11,12,16,17)])
med_corrs2_pvals
colnames(med_corrs2) <- c("AB: Attention bias threat",
                          "ER: Adaptation to emotional conflict",
                          "ER: Stroop - fear",
                          "ER: Stroop - happy",
                          "ToM: Accuracy on affective trials",
                          "ToM: Accuracy on cognitive trials",
                          "FC: Skin conductance response to CS+ vs CS-",
                          "PT: Tanner stage",
                          "IC: Reaction time on inhibit trials",
                          "IC: Reaction time on switch trials",
                          "IC: Accuracy on Stroop task",
                          "LA: Language ability",
                          "R: Reasoning",
                          "RS: Reaction time on high- vs low-reward trials",
                          "RS: Total stars")
rownames(med_corrs2) <- c("AB: Attention bias threat",
                          "ER: Adaptation to emotional conflict",
                          "ER: Stroop - fear",
                          "ER: Stroop - happy",
                          "ToM: Accuracy on affective trials",
                          "ToM: Accuracy on cognitive trials",
                          "FC: Skin conductance response to CS+ vs CS-",
                          "PT: Tanner stage",
                          "IC: Reaction time on inhibit trials",
                          "IC: Reaction time on switch trials",
                          "IC: Accuracy on Stroop task",
                          "LA: Language ability",
                          "R: Reasoning",
                          "RS: Reaction time on high- vs low-reward trials",
                          "RS: Total stars")

kbl(med_corrs2,
    caption = "Correlations between candidate mediators") %>% 
  kable_styling()

write.csv(med_corrs2, file=here("results","med_corrs_mediators_INT_EXT.csv"))

#~~ Distributions of mediators:

## Create a variable list
names(d)
vars <- names(d[,c(25:37,42,43)])
vars
med_distribs <- CreateTableOne(vars = vars,
                               strata = c("SEX"),
                               includeNA = TRUE, 
                               addOverall = TRUE,
                               data = d)
med_distribs <- data.frame(print(med_distribs, missing=TRUE))
med_distribs <- tibble::rownames_to_column(med_distribs, "Characteristic")

# Column labels
names(med_distribs) <- c("Characteristic", "Overall", "Male", "Female", "p-value", "test", "% Missing")


# Row labels
table_rows <- c("n",  
                "AB: Attention bias threat",
                "ER: Adaptation to emotional conflict",
                "ER: Stroop - fear",
                "ER: Stroop - happy",
                "ToM: Accuracy on affective trials",
                "ToM: Accuracy on cognitive trials",
                "FC: Skin conductance response to CS+ vs CS-",
                "PT: Tanner stage",
                "AL: Language ability",
                "AR: Reasoning ability",
                "IC: Reaction time on inhibit trials",
                "IC: Reaction time on switch trials",
                "IC: Accuracy on Stroop task",
                "RS: Reaction time on high- vs low-reward trials",
                "RS: Total stars")
med_distribs$Characteristic <- table_rows
kbl(med_distribs,
    caption = "Distributions of candidate mediator variables") %>% 
  kable_styling()

write.csv(med_distribs, file=here("results","Distribution_mediators_by_sex_INT_EXT.csv"))
