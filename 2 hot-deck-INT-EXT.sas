libname d "C:/Users/eks399/OneDrive - Harvard University/VDI/DT/data";

proc import datafile="C:/Users/eks399/OneDrive - Harvard University/VDI/DT/data/dat_pre_imp.csv"
        out=dat_pre_imp
        dbms=csv
        replace;
run;

proc contents data=dat_pre_imp order=varnum;run;

proc means data=dat_pre_imp;
	var num_sib first_born early_divorce num_moves preschool preg_health preg_stress preg_alc
	preterm_lowbw cesd_mom_max FIL_THREAT FIL_DEPRIVATION thr_quart thr_none dep_quart
	S1AGE SEX POV_CHRONICITY INC_NEEDS DEGREEP1P2 max_problems BIO_DAD attentothreat
	ADAPTATION STROOP_FEAR STROOP_HAPPY ACC_ATOM ACC_CTOM scr TANNER_STAGE wasitscv
	wasitscm INHIBITION_INHIBIT_BASELINE_RT INHIBITION_SWITCH_BASELINE_RT STROOP_ACC 
	BothRuns_All_Go_Trials_Accuracy BothRuns_All_NoGo_Trials_Accurac BothRuns_All_Accurate_Go_Trials_ BothRuns_All_Iccurate_NoGo_Trial
	rs_rt TotalStars INT EXT;
run;

/* E/M imputation to get parameter estimates */
proc mi data=dat_pre_imp nimpute=0;
     var num_sib first_born early_divorce num_moves preschool preg_health preg_stress preg_alc
	preterm_lowbw cesd_mom_max FIL_THREAT FIL_DEPRIVATION thr_quart thr_none dep_quart
	S1AGE SEX POV_CHRONICITY INC_NEEDS DEGREEP1P2 max_problems BIO_DAD attentothreat
	ADAPTATION STROOP_FEAR STROOP_HAPPY ACC_ATOM ACC_CTOM scr TANNER_STAGE wasitscv
	wasitscm INHIBITION_INHIBIT_BASELINE_RT INHIBITION_SWITCH_BASELINE_RT STROOP_ACC 
	BothRuns_All_Go_Trials_Accuracy BothRuns_All_NoGo_Trials_Accurac BothRuns_All_Accurate_Go_Trials_ BothRuns_All_Iccurate_NoGo_Trial
	rs_rt TotalStars INT EXT;
     em outem=d.datem; *Produces means and covariance estimates for all the variables in the var statement - can be utilized for regression estimation, but does not provide obs-wise data;
run;
proc reg data=d.datem;
 model cbcl_ext_t2=FIL_THREAT FIL_DEPRIVATION;
run;quit;

/* Bootstrap to get standard errors for the coefficients of the model above */
proc surveyselect data=dat_pre_imp method=urs n=227 reps=1000 out=bootsamp outhits;
proc mi data=bootsamp nimpute=0 noprint;
 em outem=datem;
 by replicate;
proc reg data=datem outest=a noprint;
 model cbcl_ext_t2=FIL_THREAT FIL_DEPRIVATION;
 by replicate;
 proc means data=a std;
 var FIL_THREAT FIL_DEPRIVATION;
run;

proc contents data=datem order=varnum;run;

data d.iml_boot;
	set datem;
run;

proc print data=d.iml_boot (obs=20);
	var replicate s1age;
run;

/****************************************************************************/
/* Impute the data using hotdeck											*/

proc surveyimpute data=dat_pre_imp seed=1 method=hotdeck(selection=srswr); *Hot deck imputation with simple random sampling with replacement;
	id subid;
	var num_sib first_born early_divorce num_moves preschool preg_health preg_stress preg_alc
	preterm_lowbw cesd_mom_max FIL_THREAT FIL_DEPRIVATION thr_quart thr_none dep_quart
	S1AGE SEX POV_CHRONICITY INC_NEEDS DEGREEP1P2 max_problems BIO_DAD attentothreat
	ADAPTATION STROOP_FEAR STROOP_HAPPY ACC_ATOM ACC_CTOM scr TANNER_STAGE wasitscv
	wasitscm INHIBITION_INHIBIT_BASELINE_RT INHIBITION_SWITCH_BASELINE_RT STROOP_ACC 
	BothRuns_All_Go_Trials_Accuracy BothRuns_All_NoGo_Trials_Accurac BothRuns_All_Accurate_Go_Trials_ BothRuns_All_Iccurate_NoGo_Trial
	rs_rt TotalStars INT EXT;
	output out = d.dat_hotdeck_int_ext;
run; 

proc means data=d.dat_hotdeck_int_ext;
	var num_sib first_born early_divorce num_moves preschool preg_health preg_stress preg_alc
	preterm_lowbw cesd_mom_max FIL_THREAT FIL_DEPRIVATION thr_quart thr_none dep_quart
	S1AGE SEX POV_CHRONICITY INC_NEEDS DEGREEP1P2 max_problems BIO_DAD attentothreat
	ADAPTATION STROOP_FEAR STROOP_HAPPY ACC_ATOM ACC_CTOM scr TANNER_STAGE wasitscv
	wasitscm INHIBITION_INHIBIT_BASELINE_RT INHIBITION_SWITCH_BASELINE_RT STROOP_ACC 
	BothRuns_All_Go_Trials_Accuracy BothRuns_All_NoGo_Trials_Accurac BothRuns_All_Accurate_Go_Trials_ BothRuns_All_Iccurate_NoGo_Trial
	rs_rt TotalStars INT EXT;
run;
