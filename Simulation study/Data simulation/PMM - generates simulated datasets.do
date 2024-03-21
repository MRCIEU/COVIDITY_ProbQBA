/*************************************************************************************************************************************
                                  SIMULATES DATASETS ACCORDING TO PATTERN-MIXTURE MODEL FACTORISATION
						
						SIMULATES DATASETS OF ONE SAMPLE SIZE: 100,000 OBSERVATIONS
								  
DO-FILE ARGUMENTS:
	causal:
	                !=1 SIMULATES DATASETS FOR NULL SCENARIO 
			        ==1 SIMULATES DATASETS FOR NOT NULL SCENARIO
			
	numsimdatasets: NUMBER OF SIMULATED DATASETS GENERATED						
*************************************************************************************************************************************/
args causal numsimdatasets

if `causal'==1 local causaltxt `"causaleffect"'
else local causaltxt `"nulleffect"'

* SET THE TRUE VALUES OF THE DATA GENERATION MODELS
run "Do files\PMM - sets true values of data generation models.do" `causal'

foreach samplesize of numlist 100000 {
	forvalues dataset=1(1)`numsimdatasets' {
	quietly {   
		* SIMULATE A DATASET
		run "Do files\PMM - simulates a dataset.do" `samplesize'

		* DROPS UNWANTED VARIABLES 
		drop pr_*
		
		* SAVE SIMULATED DATASETS
		export delimited using "Simulated data/PMM/`samplesize'/`causaltxt'\Dataset_`dataset'.csv", replace nolabel
	
	} // END OF QUIETLY STATEMENT
	} // END OF dataset FOR-LOOP		
} // END OF samplesize FOR-LOOP

* END OF DO-FILE