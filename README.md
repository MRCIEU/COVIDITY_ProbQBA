## Accounting for bias due to outcome data missing not at random: comparison and illustration of two approaches to probabilistic bias analysis: a simulation study

This repository contains the scripts needed to perform the probabalistic bias analyses for selection bias
due to missing data reported in the associated paper. This includes scripts to generate the simulated datasets 
and then analyse these simulations, plus scripts applying these methods to UK BioBank data.

Code for the Not-At-Random Fully Conditional Specification (NARFCS) method is available in both R (.r) and
Stata (.do), while code for the Bayesian selection model method is only available in R.

Note also that these scripts were predominantly run on the University of Bristol's High-Performance Computing
suite (https://www.bristol.ac.uk/acrc/high-performance-computing/). The shell scripts to run these jobs on this
system have not been included here.

### Simulation study

#### Data simulation
SIMULATE DATASETS OF 4 TYPES: 
(I) SELECTION MODEL FACTORISATION AND NULL EXPOSURE EFFECT, 
(II) SELECTION MODEL FACTORISATION AND NOT NULL EXPOSURE EFFECT, 
(III) PATTERN-MIXTURE MODEL FACTORISATION AND NOT NULL EXPOSURE EFFECT
(IV)  PATTERN-MIXTURE MODEL FACTORISATION AND NULL EXPOSURE EFFECT

- Master do file for simulating datasets.do
- (I AND II)
	- SM - sets true values of data generation models.do
	- SM - generates simulated datasets.do
	- SM - simulates a dataset.do
- (III AND IV)
	- PMM - sets true values of data generation models.do
	- PMM - generates simulated datasets.do
	- PMM - simulates a dataset.do

GENERATES VERY LARGE DATASET TO DETERMINE "TRUE" EXPOSURE EFFECT FOR PMM DATA MODEL
- PMM - generates large dataset for true exposure effect value.do

#### Data analysis
- Stata .do scripts to generate random seeds for each of the analyses in both the SM DGM (RandomSeeds.do) and PMM DGM (RandomSeeds_PMM.do)
- Scripts to perform the complete-cases analyses, multiple imputation, inverse-probability weighting and NARFCS analyses on the simulated datasets, and then combine the simulated results together, in both R (NARFCS_Simulations.r) and Stata (NARFCS_Simulations_Stata.do)

##### Bayesian SM
- R code to perform the Bayesian selection bias quantitative bias analysis using JAGS on the simulated datasets (BayesSM_Simulations.r)
- spreadsheet files list generated random seeds (*.csv)

### UK Biobank study
 - Script to processing the UK BioBank data ready for analysis (UKBB_Processing.r)
 - Scripts to perform the complete-cases analyses, multiple imputation, inverse-probability weighting, coding all participants without COVID testing data as 'not infected' and NARFCS analyses on the UK BioBank data, in both R (NARFCS_UKB.r) and Stata (NARFCS_UKB_Stata.do)
 - NARFCS_UKB.r - R code to perform the NARFCS quantitative bias analysis, along with other methods (complete-case analysis, inverse probability weighting [stabilised and unstabilised], standard multiple imputation, and coding all participants without COVID testing data as 'not infected')
 - R code to perform the Bayesian selection bias quantitative bias analysis using JAGS on the UK BioBank data (BayesSM_UKB.r)
