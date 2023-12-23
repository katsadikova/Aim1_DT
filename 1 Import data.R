#~~ 1 Import data.R
#------------------------------------------------
# Import .sav files and create a merged data set:
# 1.	DT_CHILD_PARENT_MASTER.sav
# 2.	DT_DOTPROBE.sav
#     - Dot probe task
#     - WASI (?) but a WASI data set sent separately as well
# 3.	DT_EMOSTROOP.sav
#     - Emotional stroop task
# 4.	DT_YSR_FU.csv
# 5.	DT_CBCL_FU.csv
# 6.	dtmhfu.csv
# 7.	dt_heirp.txt
# 8.	GNG_behavior.sav
#     - Go/No-Go task
# 9.	GO1momscr.sav
# 10.	GO2momscr.sav
# 11.	GO3momscr.sav
# 12.	GO4momscr.sav
#     - Files 9-12 are Dr. Lengua's initial cohort data (kids ages 3-6)
# 13.	DT_EF.sav
#     - Stroop task
#     - NEPSY auditory attention task
#     - NEPSY circles and squares task
#     - Various measures of executive functioning
# 14.	DT_WASI.sav
#     - Check if what's stored here matches what's in DT_DOTPROBE.sav
# 15.	DT_PINATA.sav
#     - Pinata task
# 16.	DT_CHILD_PARENT_MASTER_EMGGSR_110321.sav 
#     - Fear conditioning task
# 17.	ToM_SubjectList.csv
#     - Theory of Mind task


library(tidyverse)
library(kableExtra)
library(gtsummary)
library(expss)
library(haven)
library(sjlabelled)
library(readxl)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ Read in the data

master <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/DT_CHILD_PARENT_MASTER.sav")
dot <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/DT_DOTPROBE.sav")
emstroop <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/DT_EMOSTROOP_BASELINE.sav")
write.csv(emstroop, "/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/emstroop.csv")
emstroop <- read.csv("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/emstroop.csv")
ysr_t2 <- read.csv("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/DT_YSR_FU.csv", fileEncoding = 'UTF-8-BOM')
cbcl_t2 <- read.csv("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/DT_CBCL_FU_id.csv", fileEncoding = 'UTF-8-BOM')
int_ext_t2 <- read.csv("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/dtmhfu.csv")
pfactor_t2 <- read.csv("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/dt_heirp_fu.csv")
gonogo <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/GNG_behavior.sav")
GO1momscr <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/GO1momscr.sav")
GO2momscr <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/GO2momscr.sav")
GO3momscr <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/GO3momscr.sav")
GO4momscr <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/GO4momscr.sav")
ef <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/DT_EF.sav")
# All variables in ef are already in master, so don't need it. 
wasi <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/DT_WASI.sav")
pinata <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/DT_PINATA.sav")
fear <- read_sav("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/DT_CHILD_PARENT_MASTER_EMGGSR_110321.sav")
tom <- read_xlsx("/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/ToM_SubjectList.xlsx")

names(GO4momscr)

#~~ Pull out only CS+ Skin conductance response during acquisition from the fear data set
fear <- fear[,c("ID", "ACQ1_ASCRAmp_SR0", "ACQ1_BSCRAmp_SR0")]
#~~ Pull out only accuracy on AToM&CToM from the tom file
tom <- tom[,c("ID", "ACC_ATOM", "ACC_CTOM")]



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

#~~ Pull out only relevant variables
GO <-
  GO %>%
  dplyr::select(subid,
         m4fa004, m4fa010, m4fa005, m4fa011, 
         m4fa008, m4fa012, m4fa009, m4fa013, 
         m1fa045, m4fa048, m4fa049, m1fa096, 
         m2fa095, m3fa095, m4fa095, m1fa097,
         m2fa097, m3fa097, m4fa097, m1fa099, 
         m1fa100, m1fa102, m1fa103, m1fa120, 
         m1fa126,
         #CES-D (check that it's the parent's)
         #Pull from all 4 time points, average total CES-D scores -> 
         #Totals (m1cesd, m2cesd, m3cesd, m4cesd) need to be constructed for each time point's GO data set individually above 
         m1ce001, m1ce002, m1ce003, m1ce004,
         m1ce005, m1ce006, m1ce007, m1ce008,
         m1ce009, m1ce010, m1ce011, m1ce012,
         m1ce013, m1ce014, m1ce015, m1ce016,
         m1ce017, m1ce018, m1ce019, m1ce020,
         
         m2ce001, m2ce002, m2ce003, m2ce004,
         m2ce005, m2ce006, m2ce007, m2ce008,
         m2ce009, m2ce010, m2ce011, m2ce012,
         m2ce013, m2ce014, m2ce015, m2ce016,
         m2ce017, m2ce018, m2ce019, m2ce020,
         
         m3ce001, m3ce002, m3ce003, m3ce004,
         m3ce005, m3ce006, m3ce007, m3ce008,
         m3ce009, m3ce010, m3ce011, m3ce012,
         m3ce013, m3ce014, m3ce015, m3ce016,
         m3ce017, m3ce018, m3ce019, m3ce020,
         
         m4ce001, m4ce002, m4ce003, m4ce004,
         m4ce005, m4ce006, m4ce007, m4ce008,
         m4ce009, m4ce010, m4ce011, m4ce012,
         m4ce013, m4ce014, m4ce015, m4ce016,
         m4ce017, m4ce018, m4ce019, m4ce020)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ Merge the data sets
names(int_ext_t2)
#~~ Rename the T2 outcome variables
int_ext_t2 <-
  int_ext_t2 %>% rename(ATTENTION_t2 = ATTENTION,
                        RULE_t2 = RULE, 
                        AGGRESSIVE_t2 = AGGRESSIVE, 
                        CDI_TOT_t2 = CDI_TOT, 
                        SCARED_TOT_t2 = SCARED_TOT, 
                        PTSD_SEV_COMBINED_t2 = PTSD_SEV_COMBINED)

ysr_t2 <-
  ysr_t2 %>% rename(ysr_int_t2 = Internalizing_Problems_TScore,
                    ysr_ext_t2 = Externalizing_Problems_TScore,
                    subid = studyid) %>%
  dplyr::select(subid, ysr_int_t2, ysr_ext_t2) 

cbcl_t2 <-
  cbcl_t2 %>% rename(cbcl_int_t2 = Internalizing_Problems_TScore,
                     cbcl_ext_t2 = Externalizing_Problems_TScore,
                     subid = studyid) %>%
  dplyr::select(subid, cbcl_int_t2, cbcl_ext_t2)

#~~ Loop to rename ("ID", "id", or "Subject") as "subid"
data_list <- list(master, int_ext_t2, pfactor_t2, pinata, fear, tom)
data_list_subid <- lapply(data_list, function(data) {
  data %>%
    dplyr::rename(subid=ID)
})

gonogo <- gonogo %>%
  dplyr::rename(subid=Subject)
wasi <- wasi %>%
  dplyr::rename(subid=id)
emstroop <- rename(emstroop, subid = X.subid)

#~~ Put all data frames to be merged by subid into a list
data_list_subid <- append(data_list_subid, list(dot, gonogo, wasi, emstroop, ysr_t2, cbcl_t2))
length(data_list_subid)


dat <- Reduce(function(x, y) merge(x, y, by="subid", all=TRUE), 
       data_list_subid)
length(unique(dat$subid))

summary(as.factor(dat$subid))

dat <- merge(x=dat, y=GO, by="subid", all.x=TRUE)  

dat <- 
  dat %>%
  mutate(
    attentothreat = RT_NEUTRAL_ACCURATE - RT_ANGRY_ACCURATE
  ) 
  
save(dat, file="/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/Data/dat.Rdata")
summary(dat$Total_Problems_Total_CBCL)
summary(as.factor(GO1momscr$m1cl056g))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ Create a data dictionary showing variable and value labels

label <- as.data.frame(as.matrix(get_label(dat)))
labels <- as.data.frame(as.matrix(get_labels(dat)))
values <- as.data.frame(as.matrix(get_values(dat)))

data_dictionary <- as.data.frame(cbind(label, labels, values))
cols <- c("var_label", "value_labels", "values")
names(data_dictionary) <- cols

data_dictionary$value_labels <- vapply(data_dictionary$value_labels, paste, collapse = ", ", character(1L))
data_dictionary$values <- vapply(data_dictionary$values, paste, collapse = ", ", character(1L))

# Rerun last on 12/2 - up to date
write.csv(data_dictionary, "/Users/Kat/Dropbox/PC/Documents/PHS/Research/Deprivation and Threat/data_dictionary_alt_merge_2-10-2023.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ Create a data dictionary of the original fear data set
#~~ Many of the variables in fear that are from master are better labeled in this file

# fear <- read_sav("C:/Users/Kat/Documents/PHS/Research/Deprivation and Threat/Data/DT_CHILD_PARENT_MASTER_EMGGSR_110321.sav")
# 
# label <- as.data.frame(as.matrix(get_label(fear)))
# labels <- as.data.frame(as.matrix(get_labels(fear)))
# values <- as.data.frame(as.matrix(get_values(fear)))
# 
# data_dictionary <- as.data.frame(cbind(label, labels, values))
# cols <- c("var_label", "value_labels", "values")
# names(data_dictionary) <- cols
# 
# data_dictionary$value_labels <- vapply(data_dictionary$value_labels, paste, collapse = ", ", character(1L))
# data_dictionary$values <- vapply(data_dictionary$values, paste, collapse = ", ", character(1L))
# 
# write.csv(data_dictionary, "C:/Users/Kat/Documents/PHS/Research/Deprivation and Threat/data_dictionary_fear.csv")
