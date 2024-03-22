/*********************************************************************************************************************************                        SETS THE TRUE VALUES OF THE DATA GENERATION MODELS ACCORDING TO SELECTION MODEL FACTORISATION

	(1) JOINT MODEL FOR THE COMPLETE DATA, (2) MISSING FOR OUTCOME, (3) MISSING DATA FOR OTHER VARIABLES
	
(1) JOINT DISTRIBUTION OF CONFOUNDERS (men, degree, smoker), EXPOSURE (BMI), OUTCOME (infected) AND AUXILIARY VARIABLES (asthma, diabetes, hyperten) IS FACTORISED AS FOLLOWS: 		
					  
	p(men,age,degree,smoker,BMI,infected,asthma,diabetes,hyperten)=
		p(men)p(age│men)p(degree│age,men)p(smoker│degree,age,men)
		p(BMI│smoker,degree,age,men)p(infected│BMI,smoker,degree,age,men)
		p(asthma│infected,BMI,smoker,degree,age,men)  
		p(diabetes│asthma,infected,BMI,smoker,degree,age,men)
		p(hyperten│diabetes,asthma,infected,BMI,smoker,degree,age,men)
	
(2) MISSING DATA MECHANISM FOR THE OUTCOME - SIMULATED ACCORDING TO TWO DIFFERENT MECHANISMS: (I) NONINFORMATIVE MECHANISM AND (II) INFORMATIVE MECHANISM. IN THE SIMULATION STUDY WE DID NOT REPORT RESULTS ON THE NONINFORMATIVE MECHANISM	
	
	(I) NONINFORMATIVE SELECTION WITHOUT INTERACTIONS
	p(tested_noninform|hyperten,diabetes,asthma,sdbmi)

	(II) INFORMATIVE SELECTION WITHOUT INTERACTIONS
	p(tested_inform|hyperten,diabetes,asthma,sdbmi,infected)

(3) MISSING DATA MECHANISM FOR hyperten, sdbmi, AND smoker
p(R_hyperten|sdage,men,degree,asthma,diabetes)
p(R_smoker|sdage,men,degree,asthma,diabetes)
p(R_smoker|sdage,men,degree,asthma,diabetes)

ARGUMENTS OF THE DO-FILE
	causal
	        !=1 SETS true_infected_sdbmi=0; 
			==1 SETS true_infected_sdbmi=log(3)
***********************************************************************************************************************************/
args causal

/**********************************************************************
  (1) SET TRUE VALUES OF THE JOINT DISTRIBUTION OF THE COMPLETE DATA
**********************************************************************/
* p(men)
global true_pr_male = 0.457
global true_men_cons = logit($true_pr_male)		

* p(sdage│men)
global true_sdage_cons = -.021
global true_sdage_men = .046
global true_sdage_rmse = 1

* p(degree│sdage,men)
global true_degree_cons = -.795
global true_degree_men = .14
global true_degree_sdage = -.272

* p(smoker│degree,sdage,men)
global true_smoker_cons = -2.221
global true_smoker_men = .406
global true_smoker_sdage = -.303
global true_smoker_degree = -.6

* p(sdbmi│smoker,degree,sdage,men)
global true_sdbmi_cons = .018
global true_sdbmi_men = .168
global true_sdbmi_sdage = .024
global true_sdbmi_degree = -.26
global true_sdbmi_smoker = -.123
global true_sdbmi_rmse = .987

* p(infected│sdbmi,smoker,degree,sdage,men)
global true_infected_men = .311
global true_infected_sdage = .043
global true_infected_degree = -.318
global true_infected_smoker = .213

* p(asthma│infected,sdbmi,smoker,degree,sdage,men)  
global true_asthma_cons = -1.186
global true_asthma_men = .069
global true_asthma_sdage = .131
global true_asthma_degree = -.091
global true_asthma_smoker = .413
global true_asthma_sdbmi = .188
local obs_asthma_infected = 1.058

* p(diabetes│asthma,infected,sdbmi,smoker,degree,sdage,men)
global true_diabetes_cons = -3.097
global true_diabetes_men = .617
global true_diabetes_sdage = .426
global true_diabetes_degree = -.257
global true_diabetes_smoker = .331
global true_diabetes_sdbmi = .756
global true_diabetes_asthma = .497
local obs_diabetes_infected = .668

* p(hyperten│diabetes,asthma,infected,sdbmi,smoker,degree,sdage,men)
global true_hyperten_cons = -.895
global true_hyperten_men = .349
global true_hyperten_sdage = .649
global true_hyperten_degree = -.2
global true_hyperten_smoker = .094
global true_hyperten_sdbmi = .527
global true_hyperten_asthma = .38
global true_hyperten_diabetes = 1.374
local obs_hyperten_infected = .451

* PARAMETER VALUES THAT DIFFER ACCORDING TO NOT NULL AND NULL SCENARIOS
if `causal'==1 {
	global true_infected_cons = -3.5
	global true_infected_sdbmi = log(3)
	
	foreach name in hyperten_infected diabetes_infected asthma_infected {
		global true_`name' = log(1.75*exp(`obs_`name''))
	}
}
else {
	global true_infected_cons = -3
	
	foreach name in infected_sdbmi hyperten_infected diabetes_infected asthma_infected {
		global true_`name' = 0
	}
}

/**********************************************************************
  (2) SET TRUE VALUES OF THE MISSING DATA MECHANISM FOR THE OUTCOME
**********************************************************************/
* NONINFORMATIVE SELECTION WITHOUT INTERACTIONS
* p(tested|hyperten,diabetes,asthma,sdbmi)
local obs_testni_sdbmi = .044
local obs_testni_asthma = 1.011
local obs_testni_diabetes = .549
local obs_testni_hyperten = .452

foreach name in testni_sdbmi testni_asthma testni_diabetes testni_hyperten {
	global true_`name' = log(1.75*exp(`obs_`name''))  					
}

* INFORMATIVE SELECTION WITHOUT INTERACTIONS
* p(tested|hyperten,diabetes,asthma,sdbmi,infected)
global true_itest_sdbmi = .013
global true_itest_asthma = .991
global true_itest_diabetes = .519
global true_itest_hyperten = .441
global true_itest_infected = 7.849

* PARAMETER VALUES THAT DIFFER ACCORDING TO NOT NULL AND NULL SCENARIOS
if `causal'==1 {
	global true_testni_cons = -4.75
	global true_itest_cons = -6.3
}
else {
	global true_testni_cons = -4.6
	global true_itest_cons = -6
}
 
/**************************************************************************************
  (3) SET TRUE VALUES OF THE MISSING DATA MECHANISM FOR BMI, HYPERTENSION AND SMOKER
**************************************************************************************/
* p(Rsdbmi|men, sdage, degree, asthma, diabetes)
global true_Rsdbmi_men = .233
global true_Rsdbmi_sdage = -.057
global true_Rsdbmi_degree = -.133
global true_Rsdbmi_asthma = .363
global true_Rsdbmi_diabetes = .763
global true_Rsdbmi_cons -3.2
 
* MISSING DATA MECHANISM FOR hyperten - DEPENDS ON FULLY OBSERVED VARIABLES
global true_Rhyperten_men = -.059
global true_Rhyperten_sdage = -.029
global true_Rhyperten_degree = -.19
global true_Rhyperten_asthma = .13
global true_Rhyperten_diabetes = .192
global true_Rhyperten_cons = -2.953
 
* MISSING DATA MECHANISM FOR smoker - DEPENDS ON FULLY OBSERVED VARIABLES
global true_Rsmoker_men = .072
global true_Rsmoker_sdage = .232
global true_Rsmoker_degree = -.774
global true_Rsmoker_asthma = .169
global true_Rsmoker_diabetes = .417
global true_Rsmoker_cons = -2.9

* END OF DO FILE