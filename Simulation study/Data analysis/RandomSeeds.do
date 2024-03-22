** Generating random seeds to use in the Blue Pebble scripts
** Created 4/10/2022 by Dan Major-Smith
** Stata version 17


** Set working directory
cd "FILEPATH"


*** Causal Scenario

** IPW
clear
set obs 1000
set seed 171465

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_IPW.csv", replace


** IPW2
clear
set obs 1000
set seed 469802

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_IPW2.csv", replace


** IPW3
clear
set obs 1000
set seed 915251

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_IPW3.csv", replace


** MI
clear
set obs 1000
set seed 593160

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_MI.csv", replace


** MI_noAux
clear
set obs 1000
set seed 489795

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_MI_noAux.csv", replace


** MI_noMissingCovars
clear
set obs 1000
set seed 667215

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_MI_noMissingCovars.csv", replace


** MI_noMissingCovars_noAux
clear
set obs 1000
set seed 132366

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_MI_noMissingCovars_noAux.csv", replace


** NARMICE - Informative CSP
clear
set obs 1000
set seed 290613

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_informCSP.csv", replace


** NARMICE - Weakly-informative CSP
clear
set obs 1000
set seed 285852

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_weakInformCSP.csv", replace


** NARMICE - Non-informative CSP
clear
set obs 1000
set seed 136202

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_nonInformCSP.csv", replace


** NARMICE - Informative CSP with 5 imputations
clear
set obs 1000
set seed 898112

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_informCSP_5imps.csv", replace


** NARMICE - Weakly informative CSP with 5 imputations
clear
set obs 1000
set seed 469389

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_weakInformCSP_5imps.csv", replace


** NARMICE - Non-informative CSP with 5 imputations
clear
set obs 1000
set seed 500715

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_nonInformCSP_5imps.csv", replace


** NARMICE - Informative CSP with 5000 steps
clear
set obs 1000
set seed 603242

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_informCSP_5k.csv", replace


** NARMICE - Weakly informative CSP with 5000 steps
clear
set obs 1000
set seed 863020

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_weakInformCSP_5k.csv", replace


** NARMICE - Non-informative CSP with 5000 steps
clear
set obs 1000
set seed 792666

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_nonInformCSP_5k.csv", replace



*******************************************************************************
*** Null Scenario

** IPW
clear
set obs 1000
set seed 62082

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_IPW.csv", replace


** IPW2
clear
set obs 1000
set seed 697461

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_IPW2.csv", replace


** IPW3
clear
set obs 1000
set seed 400678

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_IPW3.csv", replace


** MI
clear
set obs 1000
set seed 977814

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_MI.csv", replace



** MI_noAux
clear
set obs 1000
set seed 269326

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_MI_noAux.csv", replace


** MI_noMissingCovars
clear
set obs 1000
set seed 582242

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_MI_noMissingCovars.csv", replace


** MI_noMissingCovars_noAux
clear
set obs 1000
set seed 502564

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_MI_noMissingCovars_noAux.csv", replace


** NARMICE - Informative CSP
clear
set obs 1000
set seed 461533

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_informCSP.csv", replace


** NARMICE - Weakly-informative CSP
clear
set obs 1000
set seed 143142

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_weakInformCSP.csv", replace


** NARMICE - Non-informative CSP
clear
set obs 1000
set seed 834778

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_nonInformCSP.csv", replace


** NARMICE - Informative CSP with 5 imputations
clear
set obs 1000
set seed 945145

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_informCSP_5imps.csv", replace


** NARMICE - Weakly informative CSP with 5 imputations
clear
set obs 1000
set seed 259553

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_weakInformCSP_5imps.csv", replace


** NARMICE - Non-informative CSP with 5 imputations
clear
set obs 1000
set seed 272431

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_nonInformCSP_5imps.csv", replace


** NARMICE - Informative CSP with 5000 steps
clear
set obs 1000
set seed 229910

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_informCSP_5k.csv", replace


** NARMICE - Weakly informative CSP with 5000 steps
clear
set obs 1000
set seed 884398

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_weakInformCSP_5k.csv", replace


** NARMICE - Non-informative CSP with 5000 steps
clear
set obs 1000
set seed 6490

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_nonInformCSP_5k.csv", replace


