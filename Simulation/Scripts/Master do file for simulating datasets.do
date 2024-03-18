/***********************************************************************************************************************************
                              MASTER DO FILE TO SIMULATE DATASETS OF 4 TYPES:
								  (I)   SELECTION MODEL FACTORISATION AND NULL EXPOSURE EFFECT
								  (II)  SELECTION MODEL FACTORISATION AND NOT NULL EXPOSURE EFFECT
								  (III) PATTERN-MIXTURE MODEL FACTORISATION AND NOT NULL EXPOSURE EFFECT
								  (IV)  PATTERN-MIXTURE MODEL FACTORISATION AND NULL EXPOSURE EFFECT
												  
***********************************************************************************************************************************/
* SET WORKING DIRECTORY
	* DO FILES STORED IN "$mydirectory\Do files\"
	* SIMULATED DATASETS STORED IN "$mydirectory\Simulated data\"
cd "$mydirectory"

/****************************************************************************
                            GENERATES 1000 DATASETS 

  IN THE SIMULATION STUDY WE DID NOT REPORT RESULTS ON SAMPLE SIZE 10,000   
*****************************************************************************/   
version 16
set seed 311021
local numsimdatasets 1000

* (I) SELECTION MODEL FACTORISATION, NULL - 10,000 AND 100,000 OBS
local causal 0
run "Do files\SM - generates simulated datasets.do" `causal' `numsimdatasets'

* (II) SELECTION MODEL FACTORISATION, NOT NULL - 10,000 AND 100,000 OBS
local causal 1
run "Do files\SM - generates simulated datasets.do" `causal' `numsimdatasets'

* (III) PATTERN-MIXUTRE MODEL FACTORISATION, NOT NULL - 100,000 OBS
set seed 1735120622
local numsimdatasets 1000
run "Do files\PMM - generates simulated datasets.do" 1 `numsimdatasets'

* (IV) PATTERN-MIXUTRE MODEL FACTORISATION, NULL - 100,000 OBS
set seed 1740120622
run "Do files\PMM - generates simulated datasets.do" 0 `numsimdatasets'

* END OF DO-FILE