##### NARFCS simulation script

### This is the R script to run through BluePebble HPC (although can easily adapy to run on standard PC/laptop)

## Created by Dan Major-Smith
## R version 4.1.0

## Load in R packages (will assume packages are already installed)
library(tidyverse)
library(survey)

# Detach the NARFCS MICE package and install/load the NARMICE package (if needed)
detach(package:mice, unload = TRUE)
library(mice)

# Note: Have installed the 'NARFCS' (or 'NARMICE') package elsewhere so does not conflict with standard 'mice' package. To install the NARFCS package, first load the 'devtools' package and then run "install_github("moreno-betancur/mice")". For more details, see: https://raw.githack.com/moreno-betancur/NARFCS/master/Vignette.html


### Some arguments from BluePebble slurm script to know dataset/loop number

# Add the args command to take arguments from the sbatch script
args <- commandArgs(trailingOnly = TRUE)

# Save args/array number as numeric
i <- as.numeric(args)
print(i)


## Set the seed - Using random seeds generated in Stata for this (each simulation will have a different set of seeds)
df_seed <- read_csv("/FILEPATH/SEED_FILE.csv")
seed <- as.numeric(df_seed[i, ])
print(seed)
set.seed(seed)
rm(df_seed)


###############################################################################################
### Read in the dataset for each array/loop in BluePebble

## This will differ, depending on whether the selection model or pattern-mixture model simulations are used, and whether the causal or null datasets (For the association between BMI and COVID testing) are used
data <- read_csv(paste0("/FILEPATH/100k_causal_csv/Dataset_", i, ".csv"))
#data <- read_csv(paste0("/FILEPATH/100k_null_csv/Dataset_", i, ".csv"))
#data <- read_csv(paste0("/FILEPATH/PMM_100k_causal_csv/Dataset_", i, ".csv"))
#data <- read_csv(paste0("/FILEPATH/PMM_100k_null_csv/Dataset_", i, ".csv"))

# Check data, and convert binary variables to factors (need this for MICE to work)
data <- data %>%
  mutate(men = as.factor(men)) %>%
  mutate(degree = as.factor(degree)) %>%
  mutate(smoker = as.factor(smoker)) %>%
  mutate(infected = as.factor(infected)) %>%
  mutate(asthma = as.factor(asthma)) %>%
  mutate(diabetes = as.factor(diabetes)) %>%
  mutate(hyperten = as.factor(hyperten)) %>%
  mutate(testni = as.factor(testni)) %>%
  mutate(itest = as.factor(itest)) %>%
  mutate(Rsdbmi = as.factor(Rsdbmi)) %>%
  mutate(Rhyperten = as.factor(Rhyperten)) %>%
  mutate(incomplete_hyperten = as.factor(incomplete_hyperten)) %>%
  mutate(Rsmoker = as.factor(Rsmoker)) %>%
  mutate(incomplete_smoker = as.factor(incomplete_smoker))

# Quick check of the data
head(data)
glimpse(data)
summary(data)



###########################################################################################
### Full data analysis and store results
mod <- glm(infected ~ men + sdage + degree + smoker + sdbmi, 
           family = "binomial", data = data)

summary(mod)

df <- as.data.frame(cbind(sim = i,
                          beta = summary(mod)$coefficients["sdbmi", "Estimate"],
                          se = summary(mod)$coefficients["sdbmi", "Std. Error"],
                          or = exp(summary(mod)$coefficients["sdbmi", "Estimate"]),
                          or_lci = exp(confint(mod)["sdbmi", "2.5 %"]),
                          or_uci = exp(confint(mod)["sdbmi", "97.5 %"])))


# Now export results to CSV
write_csv(df, paste0("/FILEPATH/Full/results_", i, ".csv"))



###########################################################################################
### Complete case analysis and store results
mod <- glm(infected ~ men + sdage + degree + incomplete_smoker + incomplete_sdbmi, 
           family = "binomial", data = data, subset = itest == 1)

summary(mod)

df <- as.data.frame(cbind(sim = i,
                          beta = summary(mod)$coefficients["incomplete_sdbmi", "Estimate"],
                          se = summary(mod)$coefficients["incomplete_sdbmi", "Std. Error"],
                          or = exp(summary(mod)$coefficients["incomplete_sdbmi", "Estimate"]),
                          or_lci = exp(confint(mod)["incomplete_sdbmi", "2.5 %"]),
                          or_uci = exp(confint(mod)["incomplete_sdbmi", "97.5 %"])))


# Now export results to CSV
write_csv(df, paste0("/FILEPATH/CCA/results_", i, ".csv"))



###########################################################################################
### Perform the weighting model to get probability of being selected/tested - Selection/testing was based on BMI, hypertension, diabetes and asthma. 

## In this script, we will combine missing 'smoking' data in with the testing/selection model

## First, though, need to remove those with missing hypertension and BMI (as both need to be part of missingness model, so have to use complete cases)
data <- data %>%
  filter(complete.cases(incomplete_hyperten, incomplete_sdbmi))

# Next, Make a combined 'tested and smoking observed' binary marker
data <- data %>%
  mutate(test_smk = ifelse(itest == 1 & Rsmoker == 0, 1, 0))

## Now for the weighting model (only including causes of testing/selection here)
tested <- glm(test_smk ~ sdbmi + asthma + diabetes + hyperten,
                 family = "binomial", data = data)

# Create predicted values
data <- data %>%
  mutate(pred_test_smk = ifelse(test_smk == 1, fitted.values(tested), NA))

# Convert probabilities to inverse weights
data <- data %>%
  mutate(IPW_test_smk = 1 / pred_test_smk)


## Now run a weighted analysis model using these IP weights

## Set-up the dataset. Specify the study design (id=~1 means there are no clusters; weights specifies the inverse-probability weights)
dstrat <- svydesign(id=~1, weights=~IPW_test_smk, data=data[data$test_smk == 1, ])

# Fit the logistic model (weights are applied automatically)
mod <- svyglm(infected ~ men + sdage + degree + smoker + sdbmi, 
                      design=dstrat, family = quasibinomial) 

		  
## Store these results (and add IP weight stats to this as well)
df <- as.data.frame(cbind(sim = i,
                          beta = summary(mod)$coefficients["sdbmi", "Estimate"],
                          se = summary(mod)$coefficients["sdbmi", "Std. Error"],
                          or = exp(summary(mod)$coefficients["sdbmi", "Estimate"]),
                          or_lci = exp(confint(mod)["sdbmi", "2.5 %"]),
                          or_uci = exp(confint(mod)["sdbmi", "97.5 %"]),
						  IPW_mean = mean(data$IPW_test_smk, na.rm = TRUE),
						  IPW_sd = sd(data$IPW_test_smk, na.rm = TRUE),
						  IPW_min = min(data$IPW_test_smk, na.rm = TRUE), 
						  IPW_5 = quantile(data$IPW_test_smk, 0.025, na.rm = TRUE),
						  IPW_25 = quantile(data$IPW_test_smk, 0.25, na.rm = TRUE),
						  IPW_median = median(data$IPW_test_smk, na.rm = TRUE),
						  IPW_75 = quantile(data$IPW_test_smk, 0.75, na.rm = TRUE),
						  IPW_95 = quantile(data$IPW_test_smk, 0.975, na.rm = TRUE),
						  IPW_max = max(data$IPW_test_smk, na.rm = TRUE)))


# Now export results to CSV
write_csv(df, paste0("/FILEPATH/IPW_main/results_", i, ".csv"))


### Also create stabilised weights to try and reduce variability in the weights

# Stabilised weighting model (exposure only)
tested_sw <- glm(test_smk ~ sdbmi, family = "binomial", data = data)

# Create predicted values
data <- data %>%
  mutate(pred_test_smk_sw = ifelse(test_smk == 1, fitted.values(tested_sw), NA))

# Convert probabilities to inverse weights
data <- data %>%
  mutate(IPW_test_smk_sw = 1 / pred_test_smk_sw)
  
# Now create stabilised weights
data <- data %>%
  mutate(sw = IPW_test_smk / IPW_test_smk_sw)
  

## Now run a weighted analysis model using these stabilised IP weights

## Set-up the dataset. Specify the study design (id=~1 means there are no clusters; weights specifies the inverse-probability weights)
dstrat_sw <- svydesign(id=~1, weights=~sw, data=data[data$test_smk == 1, ])

# Fit the logistic model (weights are applied automatically)
mod_sw <- svyglm(infected ~ men + sdage + degree + smoker + sdbmi, 
                      design=dstrat_sw, family = quasibinomial) 


## Store these results (and add IP weight stats to this as well)
df_sw <- as.data.frame(cbind(sim = i,
                          beta = summary(mod_sw)$coefficients["sdbmi", "Estimate"],
                          se = summary(mod_sw)$coefficients["sdbmi", "Std. Error"],
                          or = exp(summary(mod_sw)$coefficients["sdbmi", "Estimate"]),
                          or_lci = exp(confint(mod_sw)["sdbmi", "2.5 %"]),
                          or_uci = exp(confint(mod_sw)["sdbmi", "97.5 %"]),
						  IPW_mean = mean(data$sw, na.rm = TRUE),
						  IPW_sd = sd(data$sw, na.rm = TRUE),
						  IPW_min = min(data$sw, na.rm = TRUE), 
						  IPW_5 = quantile(data$sw, 0.025, na.rm = TRUE),
						  IPW_25 = quantile(data$sw, 0.25, na.rm = TRUE),
						  IPW_median = median(data$sw, na.rm = TRUE),
						  IPW_75 = quantile(data$sw, 0.75, na.rm = TRUE),
						  IPW_95 = quantile(data$sw, 0.975, na.rm = TRUE),
						  IPW_max = max(data$sw, na.rm = TRUE)))


# Now export results to CSV
write_csv(df_sw, paste0("/FILEPATH/IPW_main/results_", i, "_sw.csv"))



###########################################################################################
### Perform the weighting model to get probability of being selected/tested - Selection was based on BMI, hypertension, diabetes and asthma. 

## In this script, we will ignore missing 'smoking' data

## First, though, need to remove those with missing hypertension, BMI and smoking data (hypertension as missing from missingness model, BMI as missing from both missingness and substantive analysis model, and smoking as missing from substantive analysis model)
data <- data %>%
  filter(complete.cases(incomplete_hyperten, incomplete_sdbmi, incomplete_smoker))

## Now for the weighting model
tested <- glm(itest ~ sdbmi + asthma + diabetes + hyperten,
                 family = "binomial", data = data)

# Create predicted values
data <- data %>%
  mutate(pred_test = ifelse(itest == 1, fitted.values(tested), NA))

# Convert probabilities to inverse weights
data <- data %>%
  mutate(IPW_test = 1 / pred_test)


## Now run a weighted analysis model using these IP weights

## Set-up the dataset. Specify the study design (id=~1 means there are no clusters; weights specifies the inverse-probability weights)
dstrat <- svydesign(id=~1, weights=~IPW_test, data=data[data$itest == 1, ])

# Fit the logistic model (weights are applied automatically)
mod <- svyglm(infected ~ men + sdage + degree + smoker + sdbmi, 
                      design=dstrat, family = quasibinomial) 


## Store these results (and add IP weight stats to this as well)
df <- as.data.frame(cbind(sim = i,
                          beta = summary(mod)$coefficients["sdbmi", "Estimate"],
                          se = summary(mod)$coefficients["sdbmi", "Std. Error"],
                          or = exp(summary(mod)$coefficients["sdbmi", "Estimate"]),
                          or_lci = exp(confint(mod)["sdbmi", "2.5 %"]),
                          or_uci = exp(confint(mod)["sdbmi", "97.5 %"]),
						  IPW_mean = mean(data$IPW_test, na.rm = TRUE),
						  IPW_sd = sd(data$IPW_test, na.rm = TRUE),
						  IPW_min = min(data$IPW_test, na.rm = TRUE), 
						  IPW_5 = quantile(data$IPW_test, 0.025, na.rm = TRUE),
						  IPW_25 = quantile(data$IPW_test, 0.25, na.rm = TRUE),
						  IPW_median = median(data$IPW_test, na.rm = TRUE),
						  IPW_75 = quantile(data$IPW_test, 0.75, na.rm = TRUE),
						  IPW_95 = quantile(data$IPW_test, 0.975, na.rm = TRUE),
						  IPW_max = max(data$IPW_test, na.rm = TRUE)))


# Now export results to CSV
write_csv(df, paste0("/FILEPATH/IPW_spec1/results_", i, ".csv"))


### Also create stabilised weights to try and reduce variability in the weights

# Stabilised weighting model (exposure only)
tested_sw <- glm(itest ~ sdbmi, family = "binomial", data = data)

# Create predicted values
data <- data %>%
  mutate(pred_test_sw = ifelse(itest == 1, fitted.values(tested_sw), NA))

# Convert probabilities to inverse weights
data <- data %>%
  mutate(IPW_test_sw = 1 / pred_test_sw)
  
# Now create stabilised weights
data <- data %>%
  mutate(sw = IPW_test / IPW_test_sw)
  
  
## Now run a weighted analysis model using these stabilised IP weights

## Set-up the dataset. Specify the study design (id=~1 means there are no clusters; weights specifies the inverse-probability weights)
dstrat_sw <- svydesign(id=~1, weights=~sw, data=data[data$itest == 1, ])

# Fit the logistic model (weights are applied automatically)
mod_sw <- svyglm(infected ~ men + sdage + degree + smoker + sdbmi, 
                      design=dstrat_sw, family = quasibinomial) 


## Store these results (and add IP weight stats to this as well)
df_sw <- as.data.frame(cbind(sim = i,
                          beta = summary(mod_sw)$coefficients["sdbmi", "Estimate"],
                          se = summary(mod_sw)$coefficients["sdbmi", "Std. Error"],
                          or = exp(summary(mod_sw)$coefficients["sdbmi", "Estimate"]),
                          or_lci = exp(confint(mod_sw)["sdbmi", "2.5 %"]),
                          or_uci = exp(confint(mod_sw)["sdbmi", "97.5 %"]),
						  IPW_mean = mean(data$sw, na.rm = TRUE),
						  IPW_sd = sd(data$sw, na.rm = TRUE),
						  IPW_min = min(data$sw, na.rm = TRUE), 
						  IPW_5 = quantile(data$sw, 0.025, na.rm = TRUE),
						  IPW_25 = quantile(data$sw, 0.25, na.rm = TRUE),
						  IPW_median = median(data$sw, na.rm = TRUE),
						  IPW_75 = quantile(data$sw, 0.75, na.rm = TRUE),
						  IPW_95 = quantile(data$sw, 0.975, na.rm = TRUE),
						  IPW_max = max(data$sw, na.rm = TRUE)))


# Now export results to CSV
write_csv(df_sw, paste0("/FILEPATH/IPW_spec1/results_", i, "_sw.csv"))



###########################################################################################
### Perform the weighting model to get probability of being selected/tested - Selection/testing was based on BMI, hypertension, diabetes and asthma. 

## In this script, we will combine missing 'smoking' data in with the testing/selection model and include additional predictors of missing smoking data in the missingness model (men, sdage and degree)

## First, though, need to remove those with missing hypertension and BMI (as both need to be part of missingness model, so have to use complete cases)
data <- data %>%
  filter(complete.cases(incomplete_hyperten, incomplete_sdbmi))

# Next, Make a combined 'tested and smoking observed' binary marker
data <- data %>%
  mutate(test_smk = ifelse(itest == 1 & Rsmoker == 0, 1, 0))

## Now for the weighting model (including causes of testing/selection AND causes of smoking missingness here)
tested <- glm(test_smk ~ sdbmi + asthma + diabetes + hyperten + men + sdage + degree,
                 family = "binomial", data = data)

# Create predicted values
data <- data %>%
  mutate(pred_test_smk = ifelse(test_smk == 1, fitted.values(tested), NA))

# Convert probabilities to inverse weights
data <- data %>%
  mutate(IPW_test_smk = 1 / pred_test_smk)


## Now run a weighted analysis model using these IP weights

## Set-up the dataset. Specify the study design (id=~1 means there are no clusters; weights specifies the inverse-probability weights)
dstrat <- svydesign(id=~1, weights=~IPW_test_smk, data=data[data$test_smk == 1, ])

# Fit the logistic model (weights are applied automatically)
mod <- svyglm(infected ~ men + sdage + degree + smoker + sdbmi, 
                      design=dstrat, family = quasibinomial) 


## Store these results (and add IP weight stats to this as well)
df <- as.data.frame(cbind(sim = i,
                          beta = summary(mod)$coefficients["sdbmi", "Estimate"],
                          se = summary(mod)$coefficients["sdbmi", "Std. Error"],
                          or = exp(summary(mod)$coefficients["sdbmi", "Estimate"]),
                          or_lci = exp(confint(mod)["sdbmi", "2.5 %"]),
                          or_uci = exp(confint(mod)["sdbmi", "97.5 %"]),
						  IPW_mean = mean(data$IPW_test_smk, na.rm = TRUE),
						  IPW_sd = sd(data$IPW_test_smk, na.rm = TRUE),
						  IPW_min = min(data$IPW_test_smk, na.rm = TRUE), 
						  IPW_5 = quantile(data$IPW_test_smk, 0.025, na.rm = TRUE),
						  IPW_25 = quantile(data$IPW_test_smk, 0.25, na.rm = TRUE),
						  IPW_median = median(data$IPW_test_smk, na.rm = TRUE),
						  IPW_75 = quantile(data$IPW_test_smk, 0.75, na.rm = TRUE),
						  IPW_95 = quantile(data$IPW_test_smk, 0.975, na.rm = TRUE),
						  IPW_max = max(data$IPW_test_smk, na.rm = TRUE)))


# Now export results to CSV
write_csv(df, paste0("/FILEPATH/IPW_spec2/results_", i, ".csv"))


### Also create stabilised weights to try and reduce variability in the weights

# Stabilised weighting model (exposure only)
tested_sw <- glm(test_smk ~ sdbmi + men + sdage + degree, family = "binomial", data = data)

# Create predicted values
data <- data %>%
  mutate(pred_test_smk_sw = ifelse(test_smk == 1, fitted.values(tested_sw), NA))

# Convert probabilities to inverse weights
data <- data %>%
  mutate(IPW_test_smk_sw = 1 / pred_test_smk_sw)
  
# Now create stabilised weights
data <- data %>%
  mutate(sw = IPW_test_smk / IPW_test_smk_sw)
  

## Now run a weighted analysis model using these stabilised IP weights

## Set-up the dataset. Specify the study design (id=~1 means there are no clusters; weights specifies the inverse-probability weights)
dstrat_sw <- svydesign(id=~1, weights=~sw, data=data[data$test_smk == 1, ])

# Fit the logistic model (weights are applied automatically)
mod_sw <- svyglm(infected ~ men + sdage + degree + smoker + sdbmi, 
                      design=dstrat_sw, family = quasibinomial) 


## Store these results (and add IP weight stats to this as well)
df_sw <- as.data.frame(cbind(sim = i,
                          beta = summary(mod_sw)$coefficients["sdbmi", "Estimate"],
                          se = summary(mod_sw)$coefficients["sdbmi", "Std. Error"],
                          or = exp(summary(mod_sw)$coefficients["sdbmi", "Estimate"]),
                          or_lci = exp(confint(mod_sw)["sdbmi", "2.5 %"]),
                          or_uci = exp(confint(mod_sw)["sdbmi", "97.5 %"]),
						  IPW_mean = mean(data$sw, na.rm = TRUE),
						  IPW_sd = sd(data$sw, na.rm = TRUE),
						  IPW_min = min(data$sw, na.rm = TRUE), 
						  IPW_5 = quantile(data$sw, 0.025, na.rm = TRUE),
						  IPW_25 = quantile(data$sw, 0.25, na.rm = TRUE),
						  IPW_median = median(data$sw, na.rm = TRUE),
						  IPW_75 = quantile(data$sw, 0.75, na.rm = TRUE),
						  IPW_95 = quantile(data$sw, 0.975, na.rm = TRUE),
						  IPW_max = max(data$sw, na.rm = TRUE)))


# Now export results to CSV
write_csv(df_sw, paste0("/FILEPATH/IPW_spec2/results_", i, "_sw.csv"))



###########################################################################################
### Next, we want to test an MI-only approach. Here we're imputing all data from participants with missing BMI, hypertension and smoking data, as well as infection status for those tested/selected.

## Make a variable for infection status which is missing if not tested
data$infected[data$itest == 0] <- NA

# Make a new dataset just for selection independent of outcome and only keep relevant variables
data_mi <- data[, c("infected", "men", "sdage", "degree", "incomplete_smoker", "incomplete_sdbmi",
                       "asthma", "diabetes", "incomplete_hyperten")]

## Set-up the imputation methods. Will impute BMI as normal distribution (rather than default PMM), and keep 'infected_ni', 'hyperten' and 'smoker' as logistic models
meth_all <- make.method(data_mi)
meth_all
meth_all[names(meth_all) == "incomplete_sdbmi"] <- "norm"
meth_all

## Make the prediction matrix
pred_all <- make.predictorMatrix(data_mi)
pred_all

## Run a test imputation to make sure all looks fine
test <- mice(data_mi, m = 5, method = meth_all, predictorMatrix = pred_all,
             print = TRUE, maxit = 0)

## Now run the imputation model. Will create 50 imputed datasets with a burn-in period of 10 iterations each
imp <- mice(data_mi, m = 50, method = meth_all, predictorMatrix = pred_all,
               print = FALSE, maxit = 10)

## Run substantive model in each imputed dataset and combine using Rubin's Rules
model <- with(imp, glm(infected ~ men + sdage + degree + incomplete_smoker + incomplete_sdbmi, 
                          family = "binomial"))
est <- pool(model)
(covid_mi_est <- summary(est, conf.int= TRUE))


## Store these results
df <- as.data.frame(cbind(sim = i,
                          beta = covid_mi_est[covid_mi_est$term == "incomplete_sdbmi", "estimate"],
                          se = covid_mi_est[covid_mi_est$term == "incomplete_sdbmi", "std.error"],
                          or = exp(covid_mi_est[covid_mi_est$term == "incomplete_sdbmi", "estimate"]),
                          or_lci = exp(covid_mi_est[covid_mi_est$term == "incomplete_sdbmi", "2.5 %"]),
                          or_uci = exp(covid_mi_est[covid_mi_est$term == "incomplete_sdbmi", "97.5 %"])))


# Now export results to CSV
write_csv(df, paste0("/FILEPATH/MI/results_", i, ".csv"))



###########################################################################################
### Next, want to run a probabilistic NARMICE analysis

### This involves:
# 1) Selecting a random value from the prior distribution for the CSP/sensitivity/bias parameter
# 2) Apply NARMICE with this CSP to generate a single imputed dataset 
# 3) Running the substantive analysis model on this imputed dataset and store the effect estimates and standard error for our exposure 
# 4) Incorporating random sampling error 
# 5) After the full round of iterations, computing the median and 95% percentiles as credible intervals to estimate the exposure effect


# Detach the original MICE package (if needed) and load the NARMICE package
detach(package:mice, unload = TRUE)
library(mice, lib.loc = "FILEPATH")


## Make a variable for infection status which is missing if not tested
data$infected[data$itest == 0] <- NA

# Drop the 'testni' marker variable where selection/tested is independent out outcome/infection, plus the complete smoking, BMI and hypertension variables.
data_narmice <- data %>%
  select(men, sdage, degree, asthma, diabetes, itest, Rsdbmi, incomplete_sdbmi, Rhyperten, incomplete_hyperten,
         Rsmoker, incomplete_smoker, infected)

# Recode the 'tested/selected' missingness marker so that 1 = missing (rather than 1 = tested, as currently is)
data_narmice <- data_narmice %>%
  mutate(itest = ifelse(itest == 1, 0, 1)) %>%
  rename(Rtest = itest) %>%
  mutate(Rtest = as.factor(Rtest))

# Set up a new prediction matrix for the new imputation
ini <- mice(data_narmice, maxit = 0, print = TRUE)

## Specify the prediction matrix (will include all variables here, as no admin vars and don't want to exclude any)
pred_narmice <- ini$predictorMatrix
pred_narmice

# Make sure the missingness markers don't predict the variables they represent
pred_narmice["infected", "Rtest"] <- 0
pred_narmice["incomplete_sdbmi", "Rsdbmi"] <- 0
pred_narmice["incomplete_hyperten", "Rhyperten"] <- 0
pred_narmice["incomplete_smoker", "Rsmoker"] <- 0
pred_narmice

# Set-up predictor matrix for unidentifiable part
# in this case the whole matrix is zeroes because the unidentifiable part of the imputation model contains a single constant (delta*M) rather than additional contributions from the other variables in the dataset
predSens <- ini$predictorMatrix
predSens[predSens == 1] <- 0
predSens

#Set-up list with sensitivity parameter values
pSens <- rep(list(list("")), ncol(data_narmice))
names(pSens) <- names(data_narmice)
pSens

# set up vector describing manner of imputation for each variable
meth_narmice <- ini$method
meth_narmice
meth_narmice[names(meth_narmice) == "incomplete_sdbmi"] <- "norm"
meth_narmice["infected"] <- "logregSens"
meth_narmice


# Select the number of burn-in/iterations per imputations, and the number of imputations. Here, will do 10 iterations to generate 1 imputed dataset
narmice_numimps <- 1
narmice_numiter <- 10

# Number of simulations to draw bias parameter from
sims <- 10000

# To collect the parameters of interest
res <- as.data.frame(array(dim = c(sims, 23))) # Number of sensitivity values we're going to try, plus the number of parameters we're going to store (here, is 23 [see row below])
colnames(res) <- c("sim", "iteration", "csp", "msp", "beta", "se_beta", "sampprev", "beta_sampError", 
                   "cons", "cons_se", "cons_sampError", "men", "men_se", "men_sampError", 
                   "sdage", "sdage_se", "sdage_sampError", "degree", "degree_se", "degree_sampError", 
                   "smoker", "smoker_se", "smoker_sampError")
res

## Set-up prior - Select either very informative prior (mean = -7.85, SD = 1), informative prior (mean = -7.85, SD = 2) or vague prior (mean = 0, SD = 10)
prior_mean <- -7.85
#prior_mean <- 0

prior_sd <- 1
#prior_sd <- 2
#prior_sd <- 10


## Looping over each simulation, drawing a different bias parameter each time
for (j in 1:sims) {
  print(paste("Processing simulation", j, "for causal model with prior", prior_mean))
  
  # Generate random bias parameter from distribution
  bias <- rnorm(n = 1, mean = prior_mean, sd = prior_sd)
  
  # Specify a CSP value for the prediction equation for missing infection data
  pSens[["infected"]]<-list(c(bias))
  
  # NARMICE imputation
  imp_NARMICE <- mice(data_narmice, m = narmice_numimps, method = meth_narmice, predictorMatrix = pred_narmice,
                      predictorSens=predSens, parmSens=pSens,
                      print = FALSE, maxit = narmice_numiter)
  
  # Create a dataset which can be used to determine the MSP (converting the imputed dataset to wide format)
  imp_NARMICE_wide <- mice::complete(imp_NARMICE, "broad", inc=TRUE)
  
  # Create a flag to indicate when COVID status was missing 
  imp_NARMICE_wide$incomp <- as.integer(ici(imp_NARMICE_wide$infected.0))
  
  # Derive the MSP for this imputation, given this CSP
  x <- glm(infected.1 ~ incomp, family = "binomial", data = imp_NARMICE_wide)
  imp_wide_NARMICE_difflogodds <- x$coefficients[2]
  
  # Derive the prevalence of infection for this imputation, given this CSP
  prev <- glm(infected.1 ~ 1, family = binomial(link = "identity"), data = imp_NARMICE_wide)
  wholesampprev <- prev$coefficients
  
  # Run the substantive analysis model
  mod <- glm(infected.1 ~ men.1 + sdage.1 + degree.1 + incomplete_smoker.1 + incomplete_sdbmi.1, 
             family = "binomial", data = imp_NARMICE_wide)
  
  # Store these estimates in the 'results' dataframe
  res[j,"sim"] <- i
  res[j,"iteration"] <- j
  res[j,"csp"] <- bias
  res[j,"msp"] <- imp_wide_NARMICE_difflogodds
  res[j,"beta"] <- summary(mod)$coefficients["incomplete_sdbmi.1", 1]
  res[j,"se_beta"] <- summary(mod)$coefficients["incomplete_sdbmi.1", 2]
  res[j,"sampprev"] <- wholesampprev
  
  # Incorporate sampling error into the parameter estimates using standard asymptotic approximation
  res[j,"beta_sampError"] <- rnorm(n = 1, mean = res[j,"beta"], sd = res[j,"se_beta"])
  
  # Store all the other coefficients too
  res[j,"cons"] <- summary(mod)$coefficients["(Intercept)", 1]
  res[j,"cons_se"] <- summary(mod)$coefficients["(Intercept)", 2]
  res[j,"cons_sampError"] <- rnorm(n = 1, mean = res[j,"cons"], sd = res[j,"cons_se"])
  res[j,"men"] <- summary(mod)$coefficients["men.12", 1]
  res[j,"men_se"] <- summary(mod)$coefficients["men.12", 2]
  res[j,"men_sampError"] <- rnorm(n = 1, mean = res[j,"men"], sd = res[j,"men_se"])
  res[j,"sdage"] <- summary(mod)$coefficients["sdage.1", 1]
  res[j,"sdage_se"] <- summary(mod)$coefficients["sdage.1", 2]
  res[j,"sdage_sampError"] <- rnorm(n = 1, mean = res[j,"sdage"], sd = res[j,"sdage_se"])
  res[j,"degree"] <- summary(mod)$coefficients["degree.12", 1]
  res[j,"degree_se"] <- summary(mod)$coefficients["degree.12", 2]
  res[j,"degree_sampError"] <- rnorm(n = 1, mean = res[j,"degree"], sd = res[j,"degree_se"])
  res[j,"smoker"] <- summary(mod)$coefficients["incomplete_smoker.12", 1]
  res[j,"smoker_se"] <- summary(mod)$coefficients["incomplete_smoker.12", 2]
  res[j,"smoker_sampError"] <- rnorm(n = 1, mean = res[j,"smoker"], sd = res[j,"smoker_se"])
  
  print(res[j,])
  print("")
}

# Check the results table
head(res)

## Store the results from each iteration
write_csv(res, paste0("/FILEPATH/results_allIterations_", i, ".csv"))

## Summarise and store these results
df <- as.data.frame(cbind(sim = i,
                          beta_median = median(res$beta_sampError),
                          beta_mean = mean(res$beta_sampError),
                          sd = sd(res$beta_sampError),
                          perc_2.5 = quantile(res$beta_sampError, 0.025),
                          perc_97.5 = quantile(res$beta_sampError, 0.975),
                          perc_25 = quantile(res$beta_sampError, 0.25),
                          perc_75 = quantile(res$beta_sampError, 0.75),
                          cons_median = median(res$cons_sampError),
                          cons_mean = mean(res$cons_sampError),
                          cons_sd = sd(res$cons_sampError),
                          cons_perc_2.5 = quantile(res$cons_sampError, 0.025),
                          cons_perc_97.5 = quantile(res$cons_sampError, 0.975),
                          cons_perc_25 = quantile(res$cons_sampError, 0.25),
                          cons_perc_75 = quantile(res$cons_sampError, 0.75),
                          men_median = median(res$men_sampError),
                          men_mean = mean(res$men_sampError),
                          men_sd = sd(res$men_sampError),
                          men_perc_2.5 = quantile(res$men_sampError, 0.025),
                          men_perc_97.5 = quantile(res$men_sampError, 0.975),
                          men_perc_25 = quantile(res$men_sampError, 0.25),
                          men_perc_75 = quantile(res$men_sampError, 0.75),
                          sdage_median = median(res$sdage_sampError),
                          sdage_mean = mean(res$sdage_sampError),
                          sdage_sd = sd(res$sdage_sampError),
                          sdage_perc_2.5 = quantile(res$sdage_sampError, 0.025),
                          sdage_perc_97.5 = quantile(res$sdage_sampError, 0.975),
                          sdage_perc_25 = quantile(res$sdage_sampError, 0.25),
                          sdage_perc_75 = quantile(res$sdage_sampError, 0.75),
                          degree_median = median(res$degree_sampError),
                          degree_mean = mean(res$degree_sampError),
                          degree_sd = sd(res$degree_sampError),
                          degree_perc_2.5 = quantile(res$degree_sampError, 0.025),
                          degree_perc_97.5 = quantile(res$degree_sampError, 0.975),
                          degree_perc_25 = quantile(res$degree_sampError, 0.25),
                          degree_perc_75 = quantile(res$degree_sampError, 0.75),
                          smoker_median = median(res$smoker_sampError),
                          smoker_mean = mean(res$smoker_sampError),
                          smoker_sd = sd(res$smoker_sampError),
                          smoker_perc_2.5 = quantile(res$smoker_sampError, 0.025),
                          smoker_perc_97.5 = quantile(res$smoker_sampError, 0.975),
                          smoker_perc_25 = quantile(res$smoker_sampError, 0.25),
                          smoker_perc_75 = quantile(res$smoker_sampError, 0.75)
                          ))


# Now export results to CSV
write_csv(df, paste0("/FILEPATH/results_", i, ".csv"))



#############################################################################
#############################################################################
#### After running these analyses, need to combine simulated results together
#### using the 'rsimsum' package. An example is given below. Note that this
#### needs to be performed in a separate script *after* all of the simulation
#### results have been stored.

# Load the 'rsimsum' package
library(rsimsum)


##### Choose a simulation to analyse (e.g., NARFCS selection model with non-null exposure-outcome association using a very informative prio/CSP)

### Read in and append all the results files
setwd("/FILEPATH")
files <- list.files(pattern = "\\.csv$")

print(paste0("Currently analysing results from: ", getwd()))

# Drop the "Models_appended.csv" and "simsum_results.csv", if included in file list - Plus the full NARMICE iteration output
#files <- files[files != "Models_appended.csv"]
#files <- files[files != "simsum_results.csv"]
files <- files[!(startsWith(files, "results_allIterations"))]

# Read in the first file
df <-  read_csv(files[1])
head(df)

# Reading each file and append to create one file
for (f in files[-1]){
  print(paste0("On file: ", f))
  df_temp <- read_csv(f)      # read the file
  df <- rbind(df, df_temp)    # append the current file
}
head(df)

# Order by simulation number
df <- df %>%
  arrange(sim)
head(df)

# Writing the appended file  
write_csv(df, "/FILEPATH/Models_appended.csv")


## Summary stats of simulations

# Median and 95% credible intervals
median(df$beta)
quantile(df$beta, c(0.025, 0.5, 0.975))

med <- round(quantile(df$beta, 0.5), 2)
lci <- round(quantile(df$beta, 0.025), 2)
uci <- round(quantile(df$beta, 0.975), 2)

# Histogram of sll simulations
hist(df$beta, col = "lightblue", xlim = c(0.4, 1.6), xlab = "BMI log-odds",
     main = paste0("Median = ", med, " (95% CI = ", lci, " to ", uci, ")"))
abline(v = log(3), lty = 3, lwd = 3, col = "red")


## Formal simulation analysis using 'rsimsum' package
simsum_res <- simsum(data = df, estvarname = "beta", se = "se", true = log(3))
summary(simsum_res)

# Extract bias, empirical SE, model-based SE, coverage, and bias-eliminated coverage - All with 95% monte carlo intervals
simsum_table <- data.frame("model" = "XXXX",
                           "numSims" = simsum_res$summ$est[simsum_res$summ$stat == "nsim"],
                           "bias" = simsum_res$summ$est[simsum_res$summ$stat == "bias"],
                           "bias_MCSE" = simsum_res$summ$mcse[simsum_res$summ$stat == "bias"],
                           "bias_lower_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "bias"] - 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "bias"]),
                           "bias_upper_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "bias"] + 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "bias"]),
                           "empiricalSE" = simsum_res$summ$est[simsum_res$summ$stat == "empse"],
                           "empiricalSE_MCSE" = simsum_res$summ$mcse[simsum_res$summ$stat == "empse"],
                           "empSE_lower_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "empse"] - 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "empse"]),
                           "empSE_upper_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "empse"] + 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "empse"]),
                           "modelSE" = simsum_res$summ$est[simsum_res$summ$stat == "modelse"],
                           "modelSE_MCSE" = simsum_res$summ$mcse[simsum_res$summ$stat == "modelse"],
                           "modSE_lower_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "modelse"] - 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "modelse"]),
                           "modSE_upper_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "modelse"] + 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "modelse"]),
                           "coverage" = simsum_res$summ$est[simsum_res$summ$stat == "cover"],
                           "coverage_MCSE" = simsum_res$summ$mcse[simsum_res$summ$stat == "cover"],
                           "cover_lower_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "cover"] - 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "cover"]),
                           "cover_upper_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "cover"] + 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "cover"]),
						   "biasElim_cover" = simsum_res$summ$est[simsum_res$summ$stat == "becover"],
						   "BEcover_MCSE" = simsum_res$summ$mcse[simsum_res$summ$stat == "becover"],
						   "BEcover_lower_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "becover"] - 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "becover"]),
						   "BEcover_upper_MCSE" = simsum_res$summ$est[simsum_res$summ$stat == "becover"] + 
                             (1.96 * simsum_res$summ$mcse[simsum_res$summ$stat == "becover"]),
						   "median_bias" = quantile(df$beta - log(3), 0.5),
						   "median_bias_lower" = quantile(df$beta - log(3), 0.025),
						   "median_bias_upper" = quantile(df$beta - log(3), 0.975))

simsum_table

# Save this results table
write_csv(simsum_table, file = "/FILEPATH/simsum_results.csv")

###########################################################################################
### Bayesian SM analysis

# for the ith simulated data set: i=1,..., 1000

# Load the 'jagsUI' package
library(jagsUI)

## This will differ, depending on whether the selection model or pattern-mixture model simulations are used, and whether the causal or null datasets (For the association between BMI and COVID testing) are used
data <- read_csv(paste0("/FILEPATH/100k_causal_csv/Dataset_", i, ".csv"))
#data <- read_csv(paste0("/FILEPATH/100k_null_csv/Dataset_", i, ".csv"))
#data <- read_csv(paste0("/FILEPATH/PMM_100k_causal_csv/Dataset_", i, ".csv"))
#data <- read_csv(paste0("/FILEPATH/PMM_100k_null_csv/Dataset_", i, ".csv"))

## Make a variable for infection status which is missing if not tested
data$infected[data$itest == 0] <- NA

## Set the seed - Using random seeds generated in Stata for this (each simulation will have a different set of seeds)
df_seed <- read_csv("/FILEPATH/SEED_FILE.csv")
seed <- as.numeric(df_seed[i, ])
print(seed)
set.seed(seed)
rm(df_seed)

    
    ##### initial values #####
    JAGSInits <- function() {list(beta0=rnorm(1, mean=0, sd=1), beta_bmi=rnorm(1, mean=0, sd=1), beta_men=rnorm(1, mean=0, sd=1), beta_age=rnorm(1, mean=0, sd=1), beta_degree=rnorm(1, mean=0, sd=1), beta_smoker=rnorm(1, mean=0, sd=1), bmi0=rnorm(1, mean=0, sd=1), bmi_men=rnorm(1, mean=0, sd=1), bmi_age=rnorm(1, mean=0, sd=1), bmi_degree=rnorm(1, mean=0, sd=1), bmi_smoker=rnorm(1, mean=0, sd=1), bmi_prec=rgamma(1, shape=0.01, rate=0.01), smoker0=rnorm(1, mean=0, sd=1), smoker_men=rnorm(1, mean=0, sd=1), smoker_age=rnorm(1, mean=0, sd=1), smoker_degree=rnorm(1, mean=0, sd=1), asthma0=rnorm(1, mean=0, sd=1), asthma_infected=rnorm(1, mean=0, sd=1), asthma_bmi=rnorm(1, mean=0, sd=1), asthma_men=rnorm(1, mean=0, sd=1), asthma_age=rnorm(1, mean=0, sd=1), asthma_degree=rnorm(1, mean=0, sd=1), asthma_smoker=rnorm(1, mean=0, sd=1), diabetes0=rnorm(1, mean=0, sd=1), diabetes_infected=rnorm(1, mean=0, sd=1), diabetes_bmi=rnorm(1, mean=0, sd=1), diabetes_men=rnorm(1, mean=0, sd=1), diabetes_age=rnorm(1, mean=0, sd=1), diabetes_degree=rnorm(1, mean=0, sd=1), diabetes_smoker=rnorm(1, mean=0, sd=1), diabetes_asthma=rnorm(1, mean=0, sd=1), hyperten0=rnorm(1, mean=0, sd=1), hyperten_infected=rnorm(1, mean=0, sd=1), hyperten_bmi=rnorm(1, mean=0, sd=1), hyperten_men=rnorm(1, mean=0, sd=1), hyperten_age=rnorm(1, mean=0, sd=1), hyperten_degree=rnorm(1, mean=0, sd=1), hyperten_smoker=rnorm(1, mean=0, sd=1), hyperten_asthma=rnorm(1, mean=0, sd=1), hyperten_diabetes=rnorm(1, mean=0, sd=1), psi0=rnorm(1, mean=0, sd=1), psi_infected=rnorm(1, mean=0, sd=1), psi_bmi=rnorm(1, mean=0, sd=1), psi_asthma=rnorm(1, mean=0, sd=1), psi_diabetes=rnorm(1, mean=0, sd=1), psi_hyperten=rnorm(1, mean=0, sd=1))}
        
        Params <- c("beta0", "beta_bmi", "beta_men", "beta_age", "beta_degree", "beta_smoker", "asthma0", "asthma_infected", "asthma_bmi", "asthma_men", "asthma_age", "asthma_degree", "asthma_smoker", "diabetes0", "diabetes_infected", "diabetes_bmi", "diabetes_men", "diabetes_age", "diabetes_degree", "diabetes_smoker", "diabetes_asthma", "hyperten0", "hyperten_infected", "hyperten_bmi", "hyperten_men", "hyperten_age", "hyperten_degree", "hyperten_smoker", "hyperten_asthma", "hyperten_diabetes", "psi0", "psi_infected", "psi_bmi", "psi_asthma", "psi_diabetes", "psi_hyperten")

writeLines("
model
{
    for (i in 1:N) {
        infected[i] ~ dbern(infected.p[i])
        logit(infected.p[i]) <- beta0 + beta_bmi * sdbmi[i] + 
            beta_men * men[i] + beta_age * sdage[i] + beta_degree * 
            degree[i] + beta_smoker * smoker[i]
        sdbmi[i] ~ dnorm(bmi0 + bmi_men * men[i] + bmi_age * 
            sdage[i] + bmi_degree * degree[i] + bmi_smoker * 
            smoker[i], bmi_prec)
        smoker[i] ~ dbern(smoker.p[i])
        logit(smoker.p[i]) <- smoker0 + smoker_men * men[i] + 
            smoker_age * sdage[i] + smoker_degree * degree[i]
        asthma[i] ~ dbern(asthma.p[i])
        logit(asthma.p[i]) <- asthma0 + asthma_infected * infected[i] + 
            asthma_bmi * sdbmi[i] + asthma_men * men[i] + asthma_age * 
            sdage[i] + asthma_degree * degree[i] + asthma_smoker * 
            smoker[i]
        diabetes[i] ~ dbern(diabetes.p[i])
        logit(diabetes.p[i]) <- diabetes0 + diabetes_infected * 
            infected[i] + diabetes_bmi * sdbmi[i] + diabetes_men * 
            men[i] + diabetes_age * sdage[i] + diabetes_degree * 
            degree[i] + diabetes_smoker * smoker[i] + diabetes_asthma * 
            asthma[i]
        hyperten[i] ~ dbern(hyperten.p[i])
        logit(hyperten.p[i]) <- hyperten0 + hyperten_infected * 
            infected[i] + hyperten_bmi * sdbmi[i] + hyperten_men * 
            men[i] + hyperten_age * sdage[i] + hyperten_degree * 
            degree[i] + hyperten_smoker * smoker[i] + hyperten_asthma * 
            asthma[i] + hyperten_diabetes * diabetes[i]
        tested[i] ~ dbern(tested.p[i])
        logit(tested.p[i]) <- psi0 + psi_infected * infected[i] + 
            psi_bmi * sdbmi[i] + psi_asthma * asthma[i] + psi_diabetes * 
            diabetes[i] + psi_hyperten * hyperten[i]
    }
    beta0 ~ dnorm(0.00000E+00, 0.01)
    beta_bmi ~ dnorm(0.00000E+00, 0.01)
    beta_men ~ dnorm(0.00000E+00, 0.01)
    beta_age ~ dnorm(0.00000E+00, 0.01)
    beta_degree ~ dnorm(0.00000E+00, 0.01)
    beta_smoker ~ dnorm(0.00000E+00, 0.01)
    bmi0 ~ dnorm(0.00000E+00, 0.01)
    bmi_men ~ dnorm(0.00000E+00, 0.01)
    bmi_age ~ dnorm(0.00000E+00, 0.01)
    bmi_degree ~ dnorm(0.00000E+00, 0.01)
    bmi_smoker ~ dnorm(0.00000E+00, 0.01)
    bmi_prec ~ dgamma(0.01, 0.01)
    smoker0 ~ dnorm(0.00000E+00, 0.01)
    smoker_men ~ dnorm(0.00000E+00, 0.01)
    smoker_age ~ dnorm(0.00000E+00, 0.01)
    smoker_degree ~ dnorm(0.00000E+00, 0.01)
    asthma0 ~ dnorm(0.00000E+00, 0.01)
    asthma_infected ~ dnorm(0.00000E+00, 0.01)
    asthma_bmi ~ dnorm(0.00000E+00, 0.01)
    asthma_men ~ dnorm(0.00000E+00, 0.01)
    asthma_age ~ dnorm(0.00000E+00, 0.01)
    asthma_degree ~ dnorm(0.00000E+00, 0.01)
    asthma_smoker ~ dnorm(0.00000E+00, 0.01)
    diabetes0 ~ dnorm(0.00000E+00, 0.01)
    diabetes_infected ~ dnorm(0.00000E+00, 0.01)
    diabetes_bmi ~ dnorm(0.00000E+00, 0.01)
    diabetes_men ~ dnorm(0.00000E+00, 0.01)
    diabetes_age ~ dnorm(0.00000E+00, 0.01)
    diabetes_degree ~ dnorm(0.00000E+00, 0.01)
    diabetes_smoker ~ dnorm(0.00000E+00, 0.01)
    diabetes_asthma ~ dnorm(0.00000E+00, 0.01)
    hyperten0 ~ dnorm(0.00000E+00, 0.01)
    hyperten_infected ~ dnorm(0.00000E+00, 0.01)
    hyperten_bmi ~ dnorm(0.00000E+00, 0.01)
    hyperten_men ~ dnorm(0.00000E+00, 0.01)
    hyperten_age ~ dnorm(0.00000E+00, 0.01)
    hyperten_degree ~ dnorm(0.00000E+00, 0.01)
    hyperten_smoker ~ dnorm(0.00000E+00, 0.01)
    hyperten_asthma ~ dnorm(0.00000E+00, 0.01)
    hyperten_diabetes ~ dnorm(0.00000E+00, 0.01)
    psi0 ~ dnorm(0.00000E+00, 0.01)
    psi_infected ~ dnorm(prior_mean, pow(prior_sd, -2))
    psi_bmi ~ dnorm(0.00000E+00, 0.01)
    psi_asthma ~ dnorm(0.00000E+00, 0.01)
    psi_diabetes ~ dnorm(0.00000E+00, 0.01)
    psi_hyperten ~ dnorm(0.00000E+00, 0.01)
}

", con="/FILEPATH/JAGSModel.txt")


    

## Set-up prior - Select either very informative prior (mean = -7.85, SD = 1), informative prior (mean = -7.85, SD = 2) or vague prior (mean = 0, SD = 10)
prior_mean <- 7.849
#prior_mean <- 0

prior_sd <- 1
#prior_sd <- 2
#prior_sd <- 10

    data_SM <- with(data, list("N"=nrow(data), "infected"=infected, "sdbmi"=incomplete_sdbmi, "men"=men, "sdage"=sdage, "degree"=degree, "smoker"=incomplete_smoker, "asthma"=asthma, "diabetes"=diabetes, "hyperten"=incomplete_hyperten, "tested"=itest, "prior_mean"=prior_mean, "prior_sd"=prior_sd))

# Set up number of burnin (5000) and total iteraitons (50000)
SM_numburnin <- 5000
SM_numiter <- 50000

    ##### Run JAGS with burn-ins #####
## And time how long this takes
start <- Sys.time()


    res <- jags(data=data_SM, inits=JAGSInits, parameters.to.save=Params, model.file="/FILEPATH/JAGSModel.txt", n.chains=1, n.iter=SM_numiter, n.burnin=SM_numburnin)


# Calculate time taken
end <- Sys.time()

end - start

## Summarise and store these results
res_SM <- NULL
for (param in Params)
{
res_SM <- as.data.frame(cbind(res_SM, median=as.vector(unlist(res$q50)[param]), mean=as.vector(unlist(res$mean)[param]), sd=as.vector(unlist(res$sd)[param]), prec_2.5=as.vector(unlist(res$q2.5)[param]), prec_97.5=as.vector(unlist(res$q97.5)[param]), prec_25=as.vector(unlist(res$q25)[param]), prec_75=as.vector(unlist(res$q75)[param])))
colnames(res_SM)[(-6:0+ncol(res_SM))] <- paste0(param, "_", colnames(res_SM)[(-6:0+ncol(res_SM))])
}

# Now export results to CSV
write_csv(res_SM, paste0("/FILEPATH/SMresults_", i, ".csv"))
