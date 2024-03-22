** Generating random seeds to use in the Blue Pebble scripts
** Created 4/10/2022 by Dan Major-Smith
** Stata version 17


** Set working directory
cd "FILEPATH"


*** Causal Scenario

** IPW
clear
set obs 1000
set seed 159342

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_IPW_PMM.csv", replace


** IPW2
clear
set obs 1000
set seed 430908

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_IPW2_PMM.csv", replace


** IPW3
clear
set obs 1000
set seed 793805

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_IPW3_PMM.csv", replace


** MI
clear
set obs 1000
set seed 493530

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_MI_PMM.csv", replace


** NARMICE - Informative CSP
clear
set obs 1000
set seed 933042

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_informCSP_PMM.csv", replace


** NARMICE - Weakly-informative CSP
clear
set obs 1000
set seed 579428

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_weakInformCSP_PMM.csv", replace


** NARMICE - Non-informative CSP
clear
set obs 1000
set seed 939119

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_causal_NARMICE_nonInformCSP_PMM.csv", replace



*******************************************************************************
*** Null Scenario

** IPW
clear
set obs 1000
set seed 61224

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_IPW_PMM.csv", replace


** IPW2
clear
set obs 1000
set seed 449730

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_IPW2_PMM.csv", replace


** IPW3
clear
set obs 1000
set seed 187742

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_IPW3_PMM.csv", replace


** MI
clear
set obs 1000
set seed 923752

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_MI_PMM.csv", replace


** NARMICE - Informative CSP
clear
set obs 1000
set seed 65968

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_informCSP_PMM.csv", replace


** NARMICE - Weakly-informative CSP
clear
set obs 1000
set seed 312268

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_weakInformCSP_PMM.csv", replace


** NARMICE - Non-informative CSP
clear
set obs 1000
set seed 73669

gen seed = int(runiform(1, 1000000000))
sum seed

* Check no duplicates
duplicates report

* Save seeds as CSV file
export delimited "seeds_null_NARMICE_nonInformCSP_PMM.csv", replace


