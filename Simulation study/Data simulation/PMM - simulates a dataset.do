/***********************************************************************************************************************************
                                  SIMULATES A SINGLE DATASET ACCORDING TO PATTERN-MIXTURE MODEL FACTORISATION
										  
(1) THE JOINT DISTRIBUTION IS FACTORISED AS FOLLOWS: 								  
	p(men,sdage,degree,smoker,BMI,itest, infected, asthma,diabetes,hyperten)=
		p(men)p(age│men)p(degree│age,men)p(smoker│degree,sdage,men)
		p(sdbmi│smoker,degree,sdage,men)p(itest|sdbmi,smoker,degree,sdage,men)
		p(infected│itest,sdbmi,smoker,degree,sdage,men)
		p(asthma│infected,itest,sdbmi,smoker,degree,sdage,men)  
		p(diabetes│asthma,infected,itest,sdbmi,smoker,degree,sdage,men)
		p(hyperten│diabetes,asthma,infected,itest,sdbmi,smoker,degree,sdage,men)
	
(2) MISSING DATA MECHANISMS - AS PER SELECTION MECHANISM
	p(Rsdmbi|sdage,men,degree,asthma,diabetes)
	p(Rhyperten|sdage,men,degree,asthma,diabetes)
	p(Rsmoker|sdage,men,degree,asthma,diabetes)

ARGUMENTS OF THE DO-FILE
	samplesize		NUMBER OF OBSERVATIONS IN A SINGLE SIMULATED DATASET	        
***********************************************************************************************************************************/
args samplesize

clear
set obs `samplesize'

/**********************************************************************************************************************************
 COMPLETE DATA: CONFOUNDERS (men, sdage, degree, smoker), EXPOSURE (sdbmi)								  
***********************************************************************************************************************************/
* p(men)
gen men = rbinomial(1,$true_pr_male)
label define men_lb 0 "Woman" 1 "Man", replace
label values men men_lb
label variable men "Man or woman"

* p(sdage│men)
gen sdage = $true_sdage_cons + $true_sdage_men*men + $true_sdage_rmse*rnormal()
label variable sdage "Standardized age"

* p(degree│sdage,men)
gen pr_degree = invlogit($true_degree_cons + $true_degree_sdage*sdage + $true_degree_men*men)
label variable pr_degree "Probability of having a university degree"
gen degree = rbinomial(1,pr_degree)
label define yesno_lb 0 "No" 1 "Yes", replace
label values degree yesno_lb
label variable degree "Has a university degree"

* p(smoker│degree,sdage,men)
gen pr_smoker = invlogit($true_smoker_cons + $true_smoker_sdage*sdage + $true_smoker_men*men + $true_smoker_degree*degree)
label variable pr_smoker "Probability of being a smoker"
gen smoker = rbinomial(1,pr_smoker)
label values smoker yesno_lb
label variable smoker "Is a current smoker"

* p(sdbmi│smoker,degree,sdage,men)
gen sdbmi = $true_sdbmi_cons + $true_sdbmi_sdage*sdage + $true_sdbmi_men*men + $true_sdbmi_degree*degree + $true_sdbmi_smoker*smoker + $true_sdbmi_rmse*rnormal() 
label variable sdbmi "Standardized body mass index"

/***********************************************************************************************************************************
 COMPLETE DATA: SELECTION MECHANISM (itest), OUTCOME (infected) AND AUXILIARY VARIABLES (asthma, diabetes, hyperten) 								  
***********************************************************************************************************************************/
* p(itest|sdbmi,smoker,degree,sdage,men)
gen pr_itest = invlogit($true_itest_cons + $true_itest_sdbmi*sdbmi + $true_itest_smoker*smoker + $true_itest_degree*degree + $true_itest_sdage*sdage + $true_itest_men*men)
label variable pr_itest "Informative selection: probability of being tested"
gen itest = rbinomial(1,pr_itest)
label variable itest "Informative selection: Tested or not tested" 
label values itest selection_lb

* p(infected│itest,sdbmi,smoker,degree,sdage,men)
gen pr_infected = invlogit($true_infected_cons + $true_infected_sdage*sdage + $true_infected_men*men + $true_infected_degree*degree + $true_infected_smoker*smoker + $true_infected_sdbmi*sdbmi + $true_infected_itest*itest)
label variable pr_infected "Probability of being infected"
gen infected = rbinomial(1,pr_infected)
label values infected yesno_lb
label variable infected "Was infected with SARS-CoV-2"

* p(asthma│infected,itest,sdbmi,smoker,degree,sdage,men)  
gen pr_asthma = invlogit($true_asthma_cons + $true_asthma_sdage*sdage + $true_asthma_men*men + $true_asthma_degree*degree + $true_asthma_smoker*smoker + $true_asthma_sdbmi*sdbmi + $true_asthma_itest*itest + $true_asthma_infected*infected)
label variable pr_asthma "Probability of being asthmatic"
gen asthma = rbinomial(1,pr_asthma)
label values asthma yesno_lb
label variable asthma "Diagnosied with asthma"

* p(diabetes│asthma,infected,sdbmi,smoker,degree,sdage,men)
gen pr_diabetes = invlogit($true_diabetes_cons + $true_diabetes_sdage*sdage + $true_diabetes_men*men + $true_diabetes_degree*degree + $true_diabetes_smoker*smoker + $true_diabetes_sdbmi*sdbmi + $true_diabetes_itest*itest + $true_diabetes_infected*infected + $true_diabetes_asthma*asthma)
label variable pr_diabetes "Probability of being diabetic"
gen diabetes = rbinomial(1,pr_diabetes)
label values diabetes yesno_lb
label variable diabetes "Diagnosied with diabetes"

* p(hyperten│diabetes,asthma,infected,sdbmi,smoker,degree,sdage,men)
gen pr_hyperten = invlogit($true_hyperten_cons + $true_hyperten_sdage*sdage + $true_hyperten_men*men + $true_hyperten_degree*degree + $true_hyperten_smoker*smoker + $true_hyperten_sdbmi*sdbmi + $true_hyperten_itest*itest + $true_hyperten_infected*infected + $true_hyperten_asthma*asthma + $true_hyperten_diabetes*diabetes)
label variable pr_hyperten "Probability of having hypertension"
gen hyperten = rbinomial(1,pr_hyperten)
label values hyperten yesno_lb
label variable hyperten "Diagnosied with hypertension"

/**************************************************************************************
  (2) SIMULATED MISSING DATA MECHANISM FOR BMI, HYPERTENSION AND SMOKER
**************************************************************************************/
* p(Rsdbmi|sdage,men,degree,asthma,diabetes)
gen pr_Rsdbmi = invlogit($true_Rsdbmi_cons + $true_Rsdbmi_sdage*sdage + $true_Rsdbmi_men*men + $true_Rsdbmi_degree*degree + $true_Rsdbmi_asthma*asthma + $true_Rsdbmi_diabetes*diabetes)
label variable pr_Rsdbmi "Probability of sdbmi being missing"
gen Rsdbmi = rbinomial(1,pr_Rsdbmi)
label variable Rsdbmi "Missing data indicator for sdbmi"
label define missing_lb 0 "Observed" 1 "Missing", replace
label values Rsdbmi missing_lb
gen incomplete_sdbmi = sdbmi if Rsdbmi==0
label variable incomplete_sdbmi "Incompletely observed version of sdbmi"

* p(Rhyperten|sdage,men,degree,asthma,diabetes)
gen pr_Rhyperten = invlogit($true_Rhyperten_cons + $true_Rhyperten_sdage*sdage + $true_Rhyperten_men*men + $true_Rhyperten_degree*degree + $true_Rhyperten_asthma*asthma + $true_Rhyperten_diabetes*diabetes)
label variable pr_Rhyperten "Probability of hyperten being missing"
gen Rhyperten = rbinomial(1,pr_Rhyperten)
label variable Rhyperten "Missing data indicator for hyperten"
label values Rhyperten missing_lb
gen incomplete_hyperten = hyperten if Rhyperten==0
label variable incomplete_hyperten "Incompletely observed version of hypertension"

* p(Rsmoker|sdage,men,degree,asthma,diabetes)
gen pr_Rsmoker = invlogit($true_Rsmoker_cons + $true_Rsmoker_sdage*sdage + $true_Rsmoker_men*men + $true_Rsmoker_degree*degree + $true_Rsmoker_asthma*asthma + $true_Rsmoker_diabetes*diabetes)
label variable pr_Rsmoker "Probability of smoker being missing"
gen Rsmoker = rbinomial(1,pr_Rsmoker)
label variable Rsmoker "Missing data indicator for smoker"
label values Rsmoker missing_lb
gen incomplete_smoker = smoker if Rsmoker==0
label variable incomplete_smoker "Incompletely observed version of smoker"

* END OF FILE