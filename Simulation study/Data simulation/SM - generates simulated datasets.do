/*************************************************************************************************************************************
                                  SIMULATES DATASETS ACCORDING TO SELECTION MODEL FACTORISATION
						
						SIMULATES DATASETS OF TWO SAMPLE SIZES: 10,000 OBSERVATIONS AND 100,000 OBSERVATIONS
								  
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
run "Do files\SM - sets true values of data generation models.do" `causal'

foreach samplesize of numlist 10000 100000 {
	forvalues dataset=1(1)`numsimdatasets' {
	quietly {   
		* SIMULATE A DATASET
		run "Do files\SM - simulates a dataset.do" `samplesize'

		* DROPS UNWANTED VARIABLES 
		drop pr_* testni
		
		* SAVE SIMULATED DATASETS
		export delimited using "Simulated data/`samplesize'/`causaltxt'\Dataset_`dataset'.csv", replace nolabel
	
	} // END OF QUIETLY STATEMENT
	} // END OF dataset FOR-LOOP		
} // END OF samplesize FOR-LOOP

* END OF DO-FILE