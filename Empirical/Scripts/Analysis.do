/*        
Script for paper, using UK BioBank data as real-world example - Focusing here on NARFCS/NARMICE with full UKB sample and no missing covariate data

Aim is to show how to apply a probabilistic NARMICE analysis to overcome selection bias. 

Substantive model is a logistic regression, with the form: COVID ~ BMI + age + men + smoker + degree
*/

*******************************************************************************
** 1 CCA
*******************************************************************************

* Read CSV file
import delimited "$data/Dataset.csv", clear

* Fit logistic regression model
logit infected men sdage degree smoker sdbmi

* Display model summary
estat summarize

* Create a dataframe
gen beta = _b[sdbmi]
gen se = _se[sdbmi]
gen or = exp(_b[sdbmi])
gen or_lci = exp(_b[sdbmi] - 1.96 * _se[sdbmi])
gen or_uci = exp(_b[sdbmi] + 1.96 * _se[sdbmi])

keep sim - or_uci
keep if _n==1

* Export results to CSV
export delimited "$output/results.csv", replace

*******************************************************************************

/*******************************************************************************
** IPW: Perform the weighting model to get probability of being selected/tested - Selection/testing was based on BMI, hypertension, diabetes and asthma. 

*/
*******************************************************************************/



clear
tempname memhold
	postfile `memhold' str30 method beta se lci uci double p N using "$output\IPW\ipw", replace
	

import delimited "$data\Dataset.csv", clear

** Weighting model for being tested/selected
logit tested  sdbmi sdage men degree smoker asthma diabetes hyperten
	
predict pw
	
** Now derive a variable sampwgt to assign the correct sampling weights to participants:

	gen sampwgt=1/pw
	
		logit infected sdbmi men sdage degree smoker [pw=sampwgt] if tested==1


		matrix M = r(table)
		local beta = M[1,1]
		local se = M[2,1]
		local lci = M[5,1]
		local uci = M[6,1]
		local p = M[4,1]
	
		local N=e(N)
		
		
		post `memhold' ("IPW") (`beta') (`se') (`lci') (`uci') (`p') (`N') 
	
	
		
postclose `memhold'

/****************************************************************
MI only - for informative selection
Impute all data from participants with missing BMI, hypertension and smoking data, as well as infection status for those tested/selected.
Next, we want to test an MI-only approach. Here we're imputing all data from participants with missing BMI, hypertension and smoking data, as well as infection status for those tested/selected.

*****************************************************************/

*MI for informative selection - selection does depend on outcome



import delimited "$data/Dataset.csv", clear

		
mi set flong

*First, will drop the 'tested' variable

drop tested

* Check missing data in variables with missing data
mi misstable sum sdbmi hyperten smoker


* Identify which variables in the imputation model has missing information; variables used to predict covariates should be identified as well
mi register impute sdbmi hyperten smoker 

ta _mi_miss

** Perform multiple imputation
mi impute chained ///
(regress) sdbmi ///
(logit) infected  hyperten smoker ///
= men degree sdage diabetes asthma, ///
add(5) burnin(5) rseed(1234) augment force // force will force imputation, even when independent variables have missing values
* Note: right side (after = ) must contain variables with no missing information, and therefore solely considered "predictors" of missing values
* However, all variables will be included in the model to predict each missing variable

**RUN SUBSTANTIVE MODEL
mi estimate: logit infected sdbmi men sdage degree smoker

*matrix M = r(table)
		
* Create a dataframe
gen beta = _b[sdbmi]
gen se = _se[sdbmi]
gen or = exp(_b[sdbmi])
gen or_lci = exp(_b[sdbmi] - 1.96 * _se[sdbmi])
gen or_uci = exp(_b[sdbmi] + 1.96 * _se[sdbmi])

keep beta - or_uci
keep if _n==1
	
* Export results to CSV
export delimited "$output/results.csv", replace

/*******************************************************************************
probabilistic NARMICE analysis

*# 1) Selecting a random value from the prior distribution for the CSP/sensitivity/bias parameter
# 2) Apply NARMICE with this CSP to generate a single imputed dataset 
# 3) Running the substantive analysis model on this imputed dataset and store the effect estimates and standard error for our exposure 
# 4) Incorporating random sampling error 
# 5) After the full round of iterations, computing the median and 95% percentiles as credible intervals to estimate the exposure effect

**************************************************************************/

local prior `"prior1"'
**needs to be 10000
local num_simulations 10000

* Read CSV file
import delimited "$data/Dataset.csv", clear



* MISSING DATA INDICATORS FOR INCOMPLETE_INFECTED
	* RECODE SUCH THAT M_incomplete=1 DENOTES MISSING DATA
recode tested 0=1 1=0, gen(M_incomplete)
	* CHECK
tab tested M_incomplete, m	

noisily {

tempfile temp
tempname memhold
** Keep relevant variables
keep  men sdage degree smoker sdbmi asthma diabetes hyperten M_incomplete infected rsdbmi rhyperten rsmoker

save `temp', replace

set seed 220425

** CSP values	
else if "`prior'"=="prior1" {
	local hypermean -2.6
	local hypersd 0.22
}	


	
postfile `memhold' MCsimulation csp hatbeta_sdbmi se_hatbeta_sdbmi beta_samperror using "$output/pba_`prior'_narmice_tested_UKB", replace
	
forvalues n=1(1)`num_simulations' {	
	noisily di "Processing MCsimulation `n' for UKB effect with prior `prior'"
	
	use `temp', clear
	
	
* Offset for missing data	
local csp = rnormal(`hypermean',`hypersd')
	gen msp = `csp'*M_incomplete
	
	
	* IMPUTE Y
	mi set flong
mi register regular men sdage degree asthma diabetes rsdbmi rhyperten rsmoker M_incomplete
mi register imputed infected smoker hyperten sdbmi

mi impute chained ///
	(regress, omit(rsdbmi)) sdbmi ///
	(logit, offset(msp) omit(M_incomplete)) infected ///
	(logit, omit(rhyperten)) hyperten ///
	(logit, omit(rsmoker)) smoker ///
	= men degree sdage diabetes asthma rsdbmi rhyperten rsmoker M_incomplete, ///
	add(1)  burnin(10)
	
	* COMPUTE THE BIAS-ADJUSTED ESTIMATE AND CORRESPONDING STANDARD ERROR
	keep if _mi_m==1
	logit infected  sdbmi men sdage degree smoker 
local hatbeta_sdbmi = _b[sdbmi]
local se_hatbeta_sdbmi = _se[sdbmi]
		
		
* INCORPORATE RANDOM SAMPLING ERROR USING STANDARD ASYMPTOTIC APPROXIMATION
	local beta_samperror = rnormal(`hatbeta_sdbmi', `se_hatbeta_sdbmi')
	

	post `memhold'  (`n')  (`csp')  (`hatbeta_sdbmi') (`se_hatbeta_sdbmi')  (`beta_samperror')
}
postclose `memhold'
}


****SUMMARY****
use "$output/pba_`prior'_narmice_tested_UKB", clear

tempname summ
postfile `summ' N mean sd p025 p25 p50 p75 p975 se using "$output/summary/pba_`prior'_narmice_tested_UKB_summary", replace

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
