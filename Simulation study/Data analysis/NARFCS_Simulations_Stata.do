/*******************************************************************************
*Substantive model of interest - CCA

/*
logit{Pr(infected=1)}=-5.811+ 0.311×men+0.043×sdage-0.318×degree+0.213×smoker+β_BMI×sdbmi
*/
*******************************************************************************/

// Parse command line arguments
local i : word 1 of `0'

// Convert i to numeric
local i = real("`i'", .)

// Print the argument
di "`i'"

** DATA
global data "/simulations/data/causal_100k"

**OUTPUT
global output "/simulations/output/causal_100k/CCA"

		
*******************************************************************************
** 1 CCA
*******************************************************************************

* Read CSV file
import delimited "$data/Dataset_`i'.csv", clear

* Fit logistic regression model
logit infected men sdage degree smoker sdbmi

* Display model summary
estat summarize

* Create a dataframe
gen sim = `i'
gen beta = _b[sdbmi]
gen se = _se[sdbmi]
gen or = exp(_b[sdbmi])
gen or_lci = exp(_b[sdbmi] - 1.96 * _se[sdbmi])
gen or_uci = exp(_b[sdbmi] + 1.96 * _se[sdbmi])

keep sim - or_uci
keep if _n==1

* Export results to CSV
export delimited "$output/results_`i'.csv", replace

*******************************************************************************

/*******************************************************************************
** IPW 1: Perform the weighting model to get probability of being selected/tested - Selection/testing was based on BMI, hypertension, diabetes and asthma. 

In this script, we will combine missing 'smoking' data in with the testing/selection model
*/
*******************************************************************************/



clear
tempname memhold
	postfile `memhold' str30 method dataset beta se lci uci double p N using "$output\IPW_1\ipw_fully_100k_causal_2a", replace
	
	foreach dataset of numlist 1/2 {		
	
import delimited "$data\Dataset_`dataset'.csv", clear

** First remove those with missing hypertension and BMI (as both need to be part of missingness model, so have to use complete cases)
drop if incomplete_sdbmi==. | incomplete_hyperten==.

** Next, Make a combined 'tested and smoking observed' binary marker
gen test_smk = 0
replace test_smk = 1 if itest == 1 & Rsmoker == 0

** then get the weights
logit test_smk  sdbmi asthma diabetes hyperten
	
	predict pw
	
** Now derive a variable sampwgt to assign the correct sampling weights to participants:

	gen sampwgt=1/pw
	
		logit infected sdbmi men sdage degree smoker [pw=sampwgt] if itest==1


		matrix M = r(table)
		local beta = M[1,1]
		local se = M[2,1]
		local lci = M[5,1]
		local uci = M[6,1]
		local p = M[4,1]
	
		local N=e(N)
		
		
		post `memhold' ("IPW_1")  (`dataset') (`beta') (`se') (`lci') (`uci') (`p') (`N') 
	
		}
		
postclose `memhold'


/*******************************************************************************
** IPW 2: Perform the weighting model to get probability of being selected/tested - Selection was based on BMI, hypertension, diabetes and asthma. 

In this script, we will ignore missing 'smoking' data
*/
*******************************************************************************/


clear
tempname memhold
	postfile `memhold' str30 method dataset beta se lci uci double p N using "$output\causal_100k\IPW_2\ipw_partial_100k_causal_2b", replace
	
	foreach dataset of numlist 1/2 {
	
import delimited "$data\Dataset_`dataset'.csv", clear

** First remove those with missing hypertension and BMI (as both need to be part of missingness model, so have to use complete cases)
drop if incomplete_sdbmi==. | incomplete_hyperten==. | incomplete_smoker==. 

**first get the weights
logit itest  sdbmi asthma diabetes hyperten

	
	predict pw
	
** Now derive a variable sampwgt to assign the correct sampling weights to participants:

	gen sampwgt=1/pw
	
		logit infected sdbmi  men sdage degree incomplete_smoker  [pw=sampwgt] if itest==1


		matrix M = r(table)
		local beta = M[1,1]
		local se = M[2,1]
		local lci = M[5,1]
		local uci = M[6,1]
		local p = M[4,1]
	
		local N=e(N)
		
		
		post `memhold' ("IPW_2")  (`dataset') (`beta') (`se') (`lci') (`uci') (`p') (`N') 
	
		}
postclose `memhold'

/*******************************************************************************
** IPW 3: Perform the weighting model to get probability of being selected/tested - Selection/testing was based on BMI, hypertension, diabetes and asthma. 

In this script, we will combine missing 'smoking' data in with the testing/selection model and include additional predictors of missing smoking data in the missingness model (men, sdage and degree)

*******************************************************************************/


clear
tempname memhold
	postfile `memhold' str30 method dataset beta se lci uci double p N using "$output\causal_100k\IPW_2\ipw_partial_100k_causal_2c", replace
	
	foreach dataset of numlist 1/5 {
	
import delimited "$data\Dataset_`dataset'.csv", clear


** First remove those with missing hypertension and BMI (as both need to be part of missingness model, so have to use complete cases)
drop if incomplete_sdbmi==. | incomplete_hyperten==.

** Next, Make a combined 'tested and smoking observed' binary marker
gen test_smk = 0
replace test_smk = 1 if itest == 1 & Rsmoker == 0

** then the weights
logit test_smk  sdbmi  asthma diabetes hyperten men sdage degree, or
 
				 
	predict pw
	
** Now derive a variable sampwgt to assign the correct sampling weights to participants:

	gen sampwgt=1/pw
	
		logit infected sdbmi  men sdage degree smoker  [pw=sampwgt] if test_smk ==1


		matrix M = r(table)
		local beta = M[1,1]
		local se = M[2,1]
		local lci = M[5,1]
		local uci = M[6,1]
		local p = M[4,1]
	
		local N=e(N)
		
		
		post `memhold' ("IPW_3")  (`dataset') (`beta') (`se') (`lci') (`uci') (`p') (`N') 
	

		}
postclose `memhold'

/****************************************************************
MI only - for informative selection
Impute all data from participants with missing BMI, hypertension and smoking data, as well as infection status for those tested/selected.
Next, we want to test an MI-only approach. Here we're imputing all data from participants with missing BMI, hypertension and smoking data, as well as infection status for those tested/selected.

*****************************************************************/

*MI for informative selection - selection does depend on outcome

// Parse command line arguments
local i : word 1 of `0'

// Convert i to numeric
local i = real("`i'", .)

// Print the argument
di "`i'"

** DATA
global data ""

**OUTPUT
global output ""


import delimited "$data/Dataset_`i'.csv", clear

		
mi set flong

*Make a variable for infection status which is missing if testi == 0
gen infected_i=1 if itest == 1 & infected == 1
replace infected_i=0 if itest == 1 & infected == 0

fre infected_i
* Check missing data in those variables to be used for predicting missing variables
mi misstable sum incomplete_sdbmi incomplete_hyperten incomplete_smoker


* Identify which variables in the imputation model has missing information; variables used to predict covariates should be identified as well
mi register impute incomplete_sdbmi incomplete_hyperten incomplete_smoker infected_i 

ta _mi_miss

** Perform multiple imputation
mi impute chained ///
(regress) incomplete_sdbmi ///
(logit) infected_i  incomplete_hyperten incomplete_smoker ///
= men degree sdage /*itest*/ diabetes asthma, ///
add(5) burnin(5) rseed(1234) augment force // force will force imputation, even when independent variables have missing values
* Note: right side (after = ) must contain variables with no missing information, and therefore solely considered "predictors" of missing values
* However, all variables will be included in the model to predict each missing variable

*save "impute_i_`dataset'.dta", replace

**RUN SUBSTANTIVE MODEL
mi estimate: logit infected_i incomplete_sdbmi men sdage degree incomplete_smoker

*matrix M = r(table)
		
* Create a dataframe
gen sim = `i'
gen beta = _b[incomplete_sdbmi]
gen se = _se[incomplete_sdbmi]
gen or = exp(_b[incomplete_sdbmi])
gen or_lci = exp(_b[incomplete_sdbmi] - 1.96 * _se[incomplete_sdbmi])
gen or_uci = exp(_b[incomplete_sdbmi] + 1.96 * _se[incomplete_sdbmi])

keep sim - or_uci
keep if _n==1
	
* Export results to CSV
export delimited "$output/results_`i'.csv", replace

/*******************************************************************************
probabilistic NARMICE analysis

*# 1) Selecting a random value from the prior distribution for the CSP/sensitivity/bias parameter
# 2) Apply NARMICE with this CSP to generate a single imputed dataset 
# 3) Running the substantive analysis model on this imputed dataset and store the effect estimates and standard error for our exposure 
# 4) Incorporating random sampling error 
# 5) After the full round of iterations, computing the median and 95% percentiles as credible intervals to estimate the exposure effect

**************************************************************************/
// Parse command line arguments
local i : word 1 of `0'

// Convert i to numeric
local i = real("`i'", .)

// Print the argument
di "`i'"

** DATA
global data "/simulations/data/causal_100k"

**OUTPUT
global output "/simulations/output/causal_100k/PBA/informative"

********************************************************************************
local causal 1 
local prior `"informative"'
**needs to be 10000
local num_simulations 10000

* Read CSV file
import delimited "$data/Dataset_`i'.csv", clear

local causaltxt `"causal"'


* GENERATE INCOMPLETE OUTCOME VARIABLE
gen infected_i=1 if itest == 1 & infected == 1
replace infected_i=0 if itest == 1 & infected == 0

* MISSING DATA INDICATORS FOR INCOMPLETE_INFECTED
	* RECODE SUCH THAT M_incomplete=1 DENOTES MISSING DATA
recode itest 0=1 1=0, gen(M_incomplete)
	* CHECK
tab itest M_incomplete, m	

noisily {

tempfile temp
tempname memhold
** Keep relevant variables
keep  men sdage degree incomplete_smoker incomplete_sdbmi asthma diabetes incomplete_hyperten M_incomplete infected_i rsdbmi rhyperten rsmoker

save `temp', replace

set seed 220425

** CSP values	
else if "`prior'"=="informative" {
	local hypermean -7.849
	local hypersd 1
}	




	
postfile `memhold' MCsimulation csp hatbeta_sdbmi se_hatbeta_sdbmi beta_samperror using "$output/pba_`prior'_narmice_itest_`causaltxt'_dataset`i'", replace
	
forvalues n=1(1)`num_simulations' {	
	noisily di "Processing simulation `n' for `causaltxt' effect with prior `prior'"
	
	use `temp', clear
	
	
* Offset for missing data	
local csp = rnormal(`hypermean',`hypersd')
	gen msp = `csp'*M_incomplete
	
	
	* IMPUTE Y
	mi set flong
mi register regular men sdage degree asthma diabetes rsdbmi rhyperten rsmoker M_incomplete
mi register imputed infected_i incomplete_smoker incomplete_hyperten incomplete_sdbmi

mi impute chained ///
	(regress, omit(rsdbmi)) incomplete_sdbmi ///
	(logit, offset(msp) omit(M_incomplete)) infected_i ///
	(logit, omit(rhyperten)) incomplete_hyperten ///
	(logit, omit(rsmoker)) incomplete_smoker ///
	= men degree sdage diabetes asthma rsdbmi rhyperten rsmoker M_incomplete, ///
	add(1)  burnin(10)
	
	* COMPUTE THE BIAS-ADJUSTED ESTIMATE AND CORRESPONDING STANDARD ERROR
	keep if _mi_m==1
	logit infected_i  incomplete_sdbmi men sdage degree incomplete_smoker 
local hatbeta_sdbmi = _b[incomplete_sdbmi]
local se_hatbeta_sdbmi = _se[incomplete_sdbmi]
		
		
* INCORPORATE RANDOM SAMPLING ERROR USING STANDARD ASYMPTOTIC APPROXIMATION
	local beta_samperror = rnormal(`hatbeta_sdbmi', `se_hatbeta_sdbmi')
	

	post `memhold'  (`n')  (`csp')  (`hatbeta_sdbmi') (`se_hatbeta_sdbmi')  (`beta_samperror')
}
postclose `memhold'
}


****SUMMARY****
use "$output/pba_`prior'_narmice_itest_`causaltxt'_dataset`i'", clear

tempname summ
postfile `summ' N mean sd p025 p25 p50 p75 p975 se using "$output/summary/pba_`prior'_narmice_itest_`causaltxt'_dataset`i'_summary", replace

	sum beta_samperror, d
	local N=r(N)
	local mean=r(mean)
	local sd=r(sd)
    	local p25=r(p25)
	local p50=r(p50)
	local p75=r(p75)
	ci means beta_samperror
	local se=r(se)
	_pctile beta_samperror, p(2.5 97.5)
	local p025=r(r1)
	local p975=r(r2)


	post `summ' (`N') (`mean') (`sd') (`p025') (`p25') (`p50') (`p75') (`p975') (`se') 

postclose `summ'




