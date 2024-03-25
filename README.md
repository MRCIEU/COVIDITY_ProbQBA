## COVIDITY Project - Probabilistic bias analysis for selection bias in COVID-19 testing

This repository contains the scripts needed to perform the main analyses described in the associated paper ("Accounting for bias due to outcome data missing not at random: comparison and illustration of two approaches to probabilistic bias analysis: a simulation study"). This includes scripts to generate the simulated datasets and then analyse these simulations, plus scripts applying these methods to UK Biobank data.

Code for the Not-At-Random Fully Conditional Specification (NARFCS) method is available in both R (.r) and Stata (.do), while code for the Bayesian selection model method is only available in R and JAGS.

Note also that these scripts were predominantly run on the University of Bristol's High-Performance Computing suite (https://www.bristol.ac.uk/acrc/high-performance-computing/). The shell scripts to run these jobs on this system have not been included here.

### Simulation study

In this Github repository, the notation used to represent the variables in the simulation study are different from the notation used in the manuscript.

Notation in the Github repository | Notation in the manuscript | Description
--- | --- | --- 
infected | $Y$ | Partially observed outcome of the substantive analysis of interest
sdbmi | $X$ | Partially observed exposure of the substantive analysis of interest
men | $Z_1$ | Fully observed confounders of the outcome-exposure relationship
sdage | $Z_2$ | Fully observed continuous confounders of the outcome-exposure relationship
degree | $Z_3$ | Fully observed confounders of the outcome-exposure relationship
smoker | $W$ | Partially observed confounder of the outcome-exposure relationship
asthma | $A_1$ | Fully observed auxiliary variable
diabetes | $A_2$ | Fully observed auxiliary variable
hyperten | $D$ | Partially observed auxiliary variable
itest | $M^Y$ | Missing indicator of the outcome
Rsdbmi | $M^X$ | Missing indicator of the exposure
Rsmoker | $M^W$ | Missing indicator of the partially observed confounder
Rhyperten | $M^D$ | Missing indicator of the partially observed auxiliary variable


#### Data simulation
- "Master do file for simulating datasets.do" simulates datasets of 4 types:

	- (I) selection model (SM) factorisation and not null exposure effect,

(II) selection model factorisation and null exposure effect,

(III) pattern-mixture model (PMM) factorisation and not null exposure effect,

(IV) pattern-mixture model factorisation and null exposure effect,

- (I AND II)
	- "SM - sets true values of data generation models.do"
	- "SM - generates simulated datasets.do"
	- "SM - simulates a dataset.do"
- (III AND IV)
	- "PMM - sets true values of data generation models.do"
	- "PMM - generates simulated datasets.do"
	- "PMM - simulates a dataset.do"

- "PMM - generates large dataset for true exposure effect value.do" generates very large dataset to determine "TRUE" exposure effect for PMM data model

#### Data analysis
- Stata .do scripts to generate random seeds for each of the analyses in both the SM data generating model ("RandomSeeds.do") and PMM data generating model ("RandomSeeds_PMM.do")
- Scripts to perform the complete-cases analyses, multiple imputation, inverse probability weighting and NARFCS analyses on the simulated datasets, and then combine the simulated results together, in both R ("NARFCS_Simulations.r") and Stata ("NARFCS_Simulations_Stata.do")

##### Bayesian SM
- "BayesSM_Simulations.r" contains R code to perform the Bayesian selection model using JAGS on the simulated datasets
- Spreadsheet .csv files list generated random seeds used for each of the analyses

### UK Biobank study
 - "UKBB_Processing.r" processes the UK Biobank data ready for analysis
 - Scripts to perform the complete-cases analyses, multiple imputation, inverse probability weighting, coding all participants without COVID-19 testing data as 'not infected' and NARFCS analyses on the UK Biobank data, in both R ("NARFCS_UKB.r") and Stata ("NARFCS_UKB_Stata.do")
 - "NARFCS_UKB.r" contains R code to perform the NARFCS quantitative bias analysis, along with other methods (complete-case analysis, inverse probability weighting [stabilised and unstabilised], standard multiple imputation, and coding all participants without COVID-19 testing data as 'not infected')
 - "BayesSM_UKB.r" contains R code to perform the Bayesian selection model using JAGS on the UK Biobank data
