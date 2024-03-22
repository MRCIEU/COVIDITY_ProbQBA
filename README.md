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
- (I) SELECTION MODEL FACTORISATION AND NULL EXPOSURE EFFECT, 
- (II) SELECTION MODEL FACTORISATION AND NOT NULL EXPOSURE EFFECT, 
- (III) PATTERN-MIXTURE MODEL FACTORISATION AND NOT NULL EXPOSURE EFFECT
- (IV)  PATTERN-MIXTURE MODEL FACTORISATION AND NULL EXPOSURE EFFECT
	"Do files\Master do file for simulating datasets.do"
		# (I AND II)
		"Do files\SM - sets true values of data generation models.do"
		"Do files\SM - generates simulated datasets.do"
		"Do files\SM - simulates a dataset.do"
		# (III AND IV)
		"Do files\PMM - sets true values of data generation models.do"
		"Do files\PMM - generates simulated datasets.do"
		"Do files\PMM - simulates a dataset.do"

GENERATES VERY LARGE DATASET TO DETERMINE "TRUE" EXPOSURE EFFECT FOR PMM DATA MODEL
- "Do files\PMM - generates large dataset for true exposure effect value.do"

### Analyses and results processing
 - NARFCS_Simulations.r - R code to perform the NARFCS quantitative bias analysis, along with other methods 
 (analysis on the full dataset [i.e., no missing data], complete-case analysis, inverse probability weighting 
 [stabilised and unstabilised], and standard multiple imputation)
 - NARFCS_Simulations_Stata.do - Stata code to perform the NARFCS quantitative bias analysis, along with the
 other methods described above
 - BayesSM_Simulations.r - R code to perform the Bayesian selection bias quantitative bias analysis using JAGS


### UK Biobank study


COVIDITY Project - Probabalistic bias analysis for selection bias in COVID-19 testing

#### Data cleaning

#### Analyses
 - NARFCS_UKB.r - R code to perform the NARFCS quantitative bias analysis, along with other methods 
 (complete-case analysis, inverse probability weighting [stabilised and unstabilised], standard multiple 
 imputation, and coding all participants without COVID testing data as 'not infected')
 - NARFCS_UKB_Stata.do - Stata code to perform the NARFCS quantitative bias analysis, along with the
 other methods described above
 - BayesSM_UKB.r - R code to perform the Bayesian selection bias quantitative bias analysis using JAGS
