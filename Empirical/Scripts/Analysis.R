### Script for paper, using UK BioBank data as real-world example - Focusing here on NARFCS/NARMICE with full UKB sample and no missing covariate data

## Created by Dan Major-Smith
## R version 4.1.0


### Aim is to show how to apply a probabilistic NARMICE analysis to overcome selection bias. 
### Substantive model is a logistic regression, with the form: COVID ~ BMI + age + men + smoker + degree


## Clear workspace, set working directory and load packages
rm(list = ls())

setwd("FILEPATH")

#install.packages("tidyverse")
library(tidyverse)

#install.packages("survey")
library(survey)

# Detach the NARFCS MICE package and install/load the NARMICE package (if needed)
detach(package:mice, unload = TRUE)

#install.packages("mice")
library(mice)



## Read in the full dataset with no missing covariates
dat_raw <- read_csv("UKB_FullSample_MissingOutcomeOnly_zCorrected.csv")

# Make a copy, to save having to read in original again
dat <- dat_raw


## Quick summary of data
summary(dat)


### Summary stats

## Number tested
table(dat$tested)
round((table(dat$tested) / sum(table(dat$tested))) * 100, 2)

## Number positive for SARS-CoV-2
table(dat$infected)
round((table(dat$infected) / sum(table(dat$infected))) * 100, 2)

## Proportion of different covariates

# Men
table(dat$men)
round((table(dat$men) / sum(table(dat$men))) * 100, 2)

# Smoker
table(dat$smoker)
round((table(dat$smoker) / sum(table(dat$smoker))) * 100, 2)

# Degree
table(dat$degree)
round((table(dat$degree) / sum(table(dat$degree))) * 100, 2)

# Hypertension
table(dat$hyperten)
round((table(dat$hyperten) / sum(table(dat$hyperten))) * 100, 2)

# Diabetes
table(dat$diabetes)
round((table(dat$diabetes) / sum(table(dat$diabetes))) * 100, 2)

# Asthma
table(dat$asthma)
round((table(dat$asthma) / sum(table(dat$asthma))) * 100, 2)



##########################################################################################################
###### Actual analyses


######################
### Start with CCA ###
mod_cca <- glm(infected ~ sdbmi + sdage + men + degree + smoker, family = "binomial", data = dat)
summary(mod_cca)

exp(coef(summary(mod_cca)))
exp(confint(mod_cca))

# Store these results
res_cca <- as.data.frame(cbind(mod = "CCA",
                               beta = round(summary(mod_cca)$coefficients["sdbmi", "Estimate"], 3),
                               se = round(summary(mod_cca)$coefficients["sdbmi", "Std. Error"], 3),
                               or = round(exp(summary(mod_cca)$coefficients["sdbmi", "Estimate"]), 3),
                               or_lci = round(exp(confint(mod_cca)["sdbmi", "2.5 %"]), 3),
                               or_uci = round(exp(confint(mod_cca)["sdbmi", "97.5 %"]), 3)))
res_cca

# Now export results to CSV
write_csv(res_cca, ".\\Results\\UKBB_Full_CCA_Results.csv")



###################
### Next to IPW ###

set.seed(830939)

## Weighting model for being tested/selected
mod_test <- glm(tested ~ sdbmi + sdage + men + degree + smoker + asthma + diabetes + hyperten,
                family = "binomial", data = dat)
summary(mod_test)

# Create predicted probabilities of being tested
dat <- dat %>%
  mutate(pred_test = ifelse(tested == 1, fitted.values(mod_test), NA))
summary(dat$pred_test)
hist(dat$pred_test)

# Convert these probabilities to inverse weights
dat <- dat %>%
  mutate(IPW_test = 1 / pred_test)
summary(dat$IPW_test)
hist(dat$IPW_test)

## Now run a weighted analysis model using these IP weights 

## Set-up the dataset. Specify the study design (id=~1 means there are no clusters; weights=~IPW_test specifies the inverse-probability weights)
dstrat <- svydesign(id=~1, weights=~IPW_test, data=dat[!is.na(dat$IPW_test), ])

# Fit the logistic model (weights are applied automatically)
mod_ipw <- svyglm(infected ~ sdbmi + men + sdage + degree + smoker, 
                  design=dstrat, family = quasibinomial)

summary(mod_ipw)

exp(coef(summary(mod_ipw)))
exp(confint(mod_ipw))

# Store these results
res_ipw <- as.data.frame(cbind(mod = "IPW",
                               beta = round(summary(mod_ipw)$coefficients["sdbmi", "Estimate"], 3),
                               se = round(summary(mod_ipw)$coefficients["sdbmi", "Std. Error"], 3),
                               or = round(exp(summary(mod_ipw)$coefficients["sdbmi", "Estimate"]), 3),
                               or_lci = round(exp(confint(mod_ipw)["sdbmi", "2.5 %"]), 3),
                               or_uci = round(exp(confint(mod_ipw)["sdbmi", "97.5 %"]), 3)))
res_ipw

# Now export results to CSV
write_csv(res_ipw, ".\\Results\\UKBB_Full_IPW_Results.csv")

# Drop the IPW variables
dat$pred_test <- NULL
dat$IPW_test <- NULL



#### Also perform IPW using stabilised weights

## Weighting model for being tested/selected
mod_test <- glm(tested ~ sdbmi + sdage + men + degree + smoker + asthma + diabetes + hyperten,
                family = "binomial", data = dat)
summary(mod_test)

# Create predicted probabilities of being tested
dat <- dat %>%
  mutate(pred_test = ifelse(tested == 1, fitted.values(mod_test), NA))
summary(dat$pred_test)
hist(dat$pred_test)

# Convert these probabilities to inverse weights
dat <- dat %>%
  mutate(IPW_test = 1 / pred_test)
summary(dat$IPW_test)
hist(dat$IPW_test)


## Now for the denominator of the stabilised weights (using exposures only)
mod_test_sw <- glm(tested ~ sdbmi + sdage + men + degree + smoker,
                   family = "binomial", data = dat)
summary(mod_test_sw)

# Create predicted probabilities of being tested
dat <- dat %>%
  mutate(pred_test_sw = ifelse(tested == 1, fitted.values(mod_test_sw), NA))
summary(dat$pred_test_sw)
hist(dat$pred_test_sw)

# Convert these probabilities to inverse weights
dat <- dat %>%
  mutate(IPW_test_sw = 1 / pred_test_sw)
summary(dat$IPW_test_sw)
hist(dat$IPW_test_sw)


## Now create stabilised weights
dat <- dat %>%
  mutate(sw = IPW_test / IPW_test_sw)
summary(dat$sw)
hist(dat$sw)


## Now run a weighted analysis model using these IP weights

## Set-up the dataset. Specify the study design (id=~1 means there are no clusters; weights=~sw specifies the inverse-probability weights)
dstrat_sw <- svydesign(id=~1, weights=~sw, data=dat[!is.na(dat$IPW_test), ])

# Fit the logistic model (weights are applied automatically)
mod_ipw_sw <- svyglm(infected ~ sdbmi + men + sdage + degree + smoker, 
                     design=dstrat_sw, family = quasibinomial)

summary(mod_ipw_sw)

exp(coef(summary(mod_ipw_sw)))
exp(confint(mod_ipw_sw))

# Store these results
res_ipw_sw <- as.data.frame(cbind(mod = "IPW (SW)",
                                  beta = round(summary(mod_ipw_sw)$coefficients["sdbmi", "Estimate"], 3),
                                  se = round(summary(mod_ipw_sw)$coefficients["sdbmi", "Std. Error"], 3),
                                  or = round(exp(summary(mod_ipw_sw)$coefficients["sdbmi", "Estimate"]), 3),
                                  or_lci = round(exp(confint(mod_ipw_sw)["sdbmi", "2.5 %"]), 3),
                                  or_uci = round(exp(confint(mod_ipw_sw)["sdbmi", "97.5 %"]), 3)))
res_ipw_sw

# Now export results to CSV
write_csv(res_ipw_sw, ".\\Results\\UKBB_Full_IPW_SW_Results.csv")

# Drop the IPW variables
dat$pred_test <- NULL
dat$IPW_test <- NULL
dat$pred_test_sw <- NULL
dat$IPW_test_sw <- NULL
dat$sw <- NULL



##################
### Now for MI ###

set.seed(863132)

# First, will drop the 'tested' variable for MI then convert the remaining binary variables to factors (needed for MICE to work)
dat_mi <- dat %>%
  select(-tested) %>%
  mutate(men = as.factor(men)) %>%
  mutate(degree = as.factor(degree)) %>%
  mutate(smoker = as.factor(smoker)) %>%
  mutate(infected = as.factor(infected)) %>%
  mutate(asthma = as.factor(asthma)) %>%
  mutate(diabetes = as.factor(diabetes)) %>%
  mutate(hyperten = as.factor(hyperten))

summary(dat_mi)

## Set-up the imputation methods. Will impute 'infected' as logistic model
meth_all <- make.method(dat_mi)
meth_all

## Make the prediction matrix
pred_all <- make.predictorMatrix(dat_mi)
pred_all

## Run a test imputation to make sure all looks fine and no errors
test_imp <- mice(dat_mi, m = 5, method = meth_all, predictorMatrix = pred_all,
                 print = TRUE, maxit = 0)
table(test_imp$nmis)

## Now run the imputation model. Will create 50 imputed datasets with a burn-in period of 1 iterations each (as just 1 variable to impute)
imp <- mice(dat_mi, m = 50, method = meth_all, predictorMatrix = pred_all,
            print = TRUE, maxit = 1)


## Run substantive model in each imputed dataset and combine using Rubin's Rules
model_mi <- with(imp, glm(infected ~ sdbmi + men + sdage + degree + smoker, 
                          family = "binomial"))
est_mi <- pool(model_mi)
(covid_mi_est <- summary(est_mi, conf.int= TRUE))


## Prevalence of infection over imputed datasets - Around 28%
prev <- pool(with(imp, glm(infected ~ 1, family = binomial(link = "identity"))))
summary(prev)


## Store these results
res_mi <- as.data.frame(cbind(mod = "MI",
                              beta = round(covid_mi_est[covid_mi_est$term == "sdbmi", "estimate"], 3),
                              se = round(covid_mi_est[covid_mi_est$term == "sdbmi", "std.error"], 3),
                              or = round(exp(covid_mi_est[covid_mi_est$term == "sdbmi", "estimate"]), 3),
                              or_lci = round(exp(covid_mi_est[covid_mi_est$term == "sdbmi", "2.5 %"]), 3),
                              or_uci = round(exp(covid_mi_est[covid_mi_est$term == "sdbmi", "97.5 %"]), 3)))
res_mi

# Now export results to CSV
write_csv(res_mi, ".\\Results\\UKBB_Full_MI_Results.csv")



############################################
### All non-tested coded as not infected ###

dat_missingNotInfected <- dat %>%
  mutate(infected = ifelse(tested == 0, 0, infected))

## Approx 0.3% infected
summary(dat_missingNotInfected)

mod_missingNotInfected <- glm(infected ~ sdbmi + sdage + men + degree + smoker, 
                              family = "binomial", data = dat_missingNotInfected)
summary(mod_missingNotInfected)

exp(coef(summary(mod_missingNotInfected)))
exp(confint(mod_missingNotInfected))

# Store these results
res_missingNotInfected <- as.data.frame(cbind(mod = "missingNotInfected",
                               beta = round(summary(mod_missingNotInfected)$coefficients["sdbmi", "Estimate"], 3),
                               se = round(summary(mod_missingNotInfected)$coefficients["sdbmi", "Std. Error"], 3),
                               or = round(exp(summary(mod_missingNotInfected)$coefficients["sdbmi", "Estimate"]), 3),
                               or_lci = round(exp(confint(mod_missingNotInfected)["sdbmi", "2.5 %"]), 3),
                               or_uci = round(exp(confint(mod_missingNotInfected)["sdbmi", "97.5 %"]), 3)))
res_missingNotInfected

# Now export results to CSV
write_csv(res_missingNotInfected, ".\\Results\\UKBB_Full_missingNotInfected_Results.csv")




######################################
### Probabilistic NARMICE analysis ###

# Detach the original MICE package (if needed) and install/load the NARMICE package
detach(package:mice, unload = TRUE)

#install.packages("devtools")
#library(devtools)
#install_github("moreno-betancur/mice",
#               lib = "FILEPATH/mice_test")
library(mice, lib.loc = "FILEPATH/mice_test")


## Convert the 'tested' marker so that 1 = missing (rather than 1 = tested/selected)
table(dat$tested, useNA = "ifany")

dat_narmice <- dat %>%
  mutate(Rtested = 1 - tested)

table(dat_narmice$Rtested, useNA = "ifany")
table(dat_narmice$tested, useNA = "ifany")

# Convert binary variables to factors (need this for MICE to work). Note that, unlike standard MI, NARMICE suggests inclusion of the binary missingness markers when applied.
dat_narmice <- dat_narmice %>%
  mutate(men = as.factor(men)) %>%
  mutate(degree = as.factor(degree)) %>%
  mutate(smoker = as.factor(smoker)) %>%
  mutate(infected = as.factor(infected)) %>%
  mutate(asthma = as.factor(asthma)) %>%
  mutate(diabetes = as.factor(diabetes)) %>%
  mutate(hyperten = as.factor(hyperten)) %>%
  mutate(Rtested = as.factor(Rtested)) %>%
  mutate(tested = as.factor(tested))

summary(dat_narmice)


### Before running the probabilistic NARMICE analysis, we first need to calibrate the CSP value based on either the MSP and/or known prevalence (here, will calibrate using estimated prevalence). Here, we are following algorithm 3 of Tompsett (2018), where we first use a wide range of MSPs/poplation parameters, then narrow this down until we reach our target value. We will repeat this for the main CSP value (prevalance of 3.2%), plus the upper and lower ranges of this CSP (estimated to be between 2.2% and 4.2%). 

# Now drop the original 'tested' variable, as not needed below
dat_narmice <- dat_narmice %>%
  select(-tested)


### Now, need to run the Tompsett algorithm to find the CSP linked to said prevalance estimates

## First iteration will use a wide range of CSPs, then narrow after, until reach the intended prevalence estimate within 4 decimal places

# Set up a prediction matrix for the imputation
ini <- mice(dat_narmice, maxit = 0, print = TRUE)

## Specify the prediction matrix (will include all variables here, as no admin vars and don't want to exclude any)
pred_narmice <- ini$predictorMatrix
pred_narmice

# Make sure the missingness markers don't predict the variables they represent
pred_narmice["infected", "Rtested"] <- 0
pred_narmice

# Set-up predictor matrix for unidentifiable part
# in this case the whole matrix is zeroes because the unidentifiable part of the imputation model contains a single constant (delta*M) rather than additional contributions from the other variables in the dataset
predSens <- ini$predictorMatrix
predSens[predSens == 1] <- 0
predSens

#Set-up list with sensitivity parameter values
pSens <- rep(list(list("")), ncol(dat_narmice))

names(pSens) <- names(dat_narmice)
pSens

# set up vector describing manner of imputation for each variable. Only imputing 'infected', which is a logistic model, but as we're adding a sensitivity parameter, we need to specify the 'logregSens' option
meth_narmice <- ini$method
meth_narmice
meth_narmice["infected"] <- "logregSens"
meth_narmice

# Set up number of imputations (5) and burn-in period/iterations (1; only need 1 as only 1 variable to impute)
narmice_numimps <- 5
narmice_numiter <- 1

# To collect the parameters of interest
res <- as.data.frame(array(dim = c(dim = length(seq.int(-4, -1, by = 0.25)), 3))) # Sensitivity parameters to store
colnames(res) <- c("csp", "msp", "sampprev")
res


## Loop over the different CSP values to test
k <- 0
for (j in seq.int(-4, -1, by = 0.25)) {
  k <- k+1
  print(paste0("Delta = ", j))

  # specify a delta value for the prediction equation for missing infection data
  pSens[["infected"]]<-list(c(j))

  # NARMICE imputation
  new.seed <- 230324
  imp_NARMICE <- mice(dat_narmice, m = narmice_numimps, method = meth_narmice, predictorMatrix = pred_narmice,
                      predictorSens=predSens, parmSens=pSens,
                      print = TRUE, seed = new.seed, maxit = narmice_numiter)

  # Estimate and store the MSP
  mspest <- round(summary(pool(with(imp_NARMICE, glm(formula = infected ~ Rtested,
                                                     family = "binomial"))))["Rtested1", "est"], 2)

  # derive the prevalence of infection in the whole population
  wholesampprev <- round(summary(pool(with(imp_NARMICE,
                                           glm(formula = infected ~ 1,
                                               family = binomial(link = "identity")))))["(Intercept)", "est"], 4)

  # Store these estimates in the 'tipping' dataframe
  res[k,"csp"] <- j
  res[k,"msp"] <- mspest
  res[k,"sampprev"] <- wholesampprev

}


## Check the results, to see where to focus on the second round of imputations to estimate the sample prevalence to 4 DPs
res


## For the mean CSP, will focus on CSPs between -2.5 and -2.6 (prevalence approx 3.2%)

# To collect the parameters of interest
res_prev <- as.data.frame(array(dim = c(dim = length(seq.int(-2.6, -2.5, by = 0.01)), 3))) # Sensitivity parameters to store
colnames(res_prev) <- c("csp", "msp", "sampprev")
res_prev

## Loop over the different CSP values to test
k <- 0
for (j in seq.int(-2.6, -2.5, by = 0.01)) {
  k <- k+1
  print(paste0("Delta = ", j))

  # specify a delta value for the prediction equation for missing infection data
  pSens[["infected"]]<-list(c(j))

  # NARMICE imputation - Also make a new seed for each loop, else starting values may be similar (I think?)
  new.seed <- 230324
  imp_NARMICE <- mice(dat_narmice, m = narmice_numimps, method = meth_narmice, predictorMatrix = pred_narmice,
                      predictorSens=predSens, parmSens=pSens,
                      print = TRUE, seed = new.seed, maxit = narmice_numiter)

  # Estimate and store the MSP
  mspest <- round(summary(pool(with(imp_NARMICE, glm(formula = infected ~ Rtested,
                                                     family = "binomial"))))["Rtested1", "est"], 2)

  # derive the prevalence of infection in the whole population
  wholesampprev <- round(summary(pool(with(imp_NARMICE,
                                           glm(formula = infected ~ 1,
                                               family = binomial(link = "identity")))))["(Intercept)", "est"], 4)

  # Store these estimates in the results dataframe
  res_prev[k,"csp"] <- j
  res_prev[k,"msp"] <- mspest
  res_prev[k,"sampprev"] <- wholesampprev

}

# See which CSP corresponds to the correct prevalence
res_prev

res_prev[which.min(abs(res_prev$sampprev - 0.0320)), ]
csp <- res_prev$csp[which.min(abs(res_prev$sampprev - 0.0320))]
csp


## For the lower prevalence CSP, will focus on CSPs between -2.95 and -3.05 (prevalence approx 2.2%)

# To collect the parameters of interest
res_prev_lower <- as.data.frame(array(dim = c(dim = length(seq.int(-3.05, -2.95, by = 0.01)), 3))) # Sensitivity parameters to store
colnames(res_prev_lower) <- c("csp", "msp", "sampprev")
res_prev_lower

## Loop over the different CSP values to test
k <- 0
for (j in seq.int(-3.05, -2.95, by = 0.01)) {
  k <- k+1
  print(paste0("Delta = ", j))

  # specify a delta value for the prediction equation for missing infection data
  pSens[["infected"]]<-list(c(j))

  # NARMICE imputation
  new.seed <- 230324
  imp_NARMICE <- mice(dat_narmice, m = narmice_numimps, method = meth_narmice, predictorMatrix = pred_narmice,
                      predictorSens=predSens, parmSens=pSens,
                      print = TRUE, seed = new.seed, maxit = narmice_numiter)

  # Estimate and store the MSP
  mspest <- round(summary(pool(with(imp_NARMICE, glm(formula = infected ~ Rtested,
                                                     family = "binomial"))))["Rtested1", "est"], 2)

  # derive the prevalence of infection in the whole population
  wholesampprev <- round(summary(pool(with(imp_NARMICE,
                                           glm(formula = infected ~ 1,
                                               family = binomial(link = "identity")))))["(Intercept)", "est"], 4)

  # Store these estimates in the results dataframe
  res_prev_lower[k,"csp"] <- j
  res_prev_lower[k,"msp"] <- mspest
  res_prev_lower[k,"sampprev"] <- wholesampprev

}

# See which CSP corresponds to the correct lower prevalence
res_prev_lower

res_prev_lower[which.min(abs(res_prev_lower$sampprev - 0.0220)), ]
csp_lower <- res_prev_lower$csp[res_prev_lower.min(abs(res_prev_lower$sampprev - 0.0220))]
csp_lower


## For the upper prevalence CSP, will focus on CSPs between -2.3 and -2.2 (prevalance approx 4.2%)

# To collect the parameters of interest
res_prev_upper <- as.data.frame(array(dim = c(dim = length(seq.int(-2.3, -2.2, by = 0.01)), 3))) # Sensitivity parameters to store
colnames(res_prev_upper) <- c("csp", "msp", "sampprev")
res_prev_upper

## Loop over the different CSP values to test
k <- 0
for (j in seq.int(-2.3, -2.2, by = 0.01)) {
  k <- k+1
  print(paste0("Delta = ", j))

  # specify a delta value for the prediction equation for missing infection data
  pSens[["infected"]]<-list(c(j))

  # NARMICE imputation
  new.seed <- 230324
  imp_NARMICE <- mice(dat_narmice, m = narmice_numimps, method = meth_narmice, predictorMatrix = pred_narmice,
                      predictorSens=predSens, parmSens=pSens,
                      print = TRUE, seed = new.seed, maxit = narmice_numiter)

  # Estimate and store the MSP
  mspest <- round(summary(pool(with(imp_NARMICE, glm(formula = infected ~ Rtested,
                                                     family = "binomial"))))["Rtested1", "est"], 2)

  # derive the prevalence of infection in the whole population
  wholesampprev <- round(summary(pool(with(imp_NARMICE,
                                           glm(formula = infected ~ 1,
                                               family = binomial(link = "identity")))))["(Intercept)", "est"], 4)

  # Store these estimates in the results dataframe
  res_prev_upper[k,"csp"] <- j
  res_prev_upper[k,"msp"] <- mspest
  res_prev_upper[k,"sampprev"] <- wholesampprev

}

# See which CSP corresponds to the correct upper prevalence
res_prev_upper

res_prev_upper[which.min(abs(res_prev_upper$sampprev - 0.0420)), ]
csp_upper <- res_prev_upper$csp[which.min(abs(res_prev_upper$sampprev - 0.0420))]
csp_upper


### Finally, using these CSP values, derive an approximate value for the SD of the CSP, such that 95% of values sampled from a normal distribution centered on the mean CSP value would fall within the upper and lower CSP values
csp; csp_lower; csp_upper

# Estimate SD for lower CSP value
csp_lower_sd <- (csp - csp_lower) / 1.96
csp_lower_sd

# Estimate SD for upper CSP value
csp_upper_sd <- (csp_upper - csp) / 1.96
csp_upper_sd

## The SDs differ quite considerably, as upper and lower CSP estimates not symmetrical.

# Decided to use the larger SD (from the lower CSP bound), so the coverage is wider (approx. 2.2 to 4.8%).


### Once have final mean and SD CSP, could comment out the algorithms above to save processing time - Also reduce CSP to 1DP, and SD to 2DPs
mean_csp <- -2.6
sd_csp <- 0.22



#### Next, run the probabilistic NARMICE analysis (using an informative prior/conditional sensitivity parameter).

### This involves:
# 1) Selecting a random value from the prior distribution for the CSP/sensitivity/bias parameter
# 2) Apply NARMICE with this CSP to generate a single imputed dataset 
# 3) Running the substantive analysis model on this imputed dataset and store the effect estimates and standard error for our exposure 
# 4) Incorporating random sampling error 
# 5) After the full round of iterations, computing the median and 95% percentiles as credible intervals to estimate the exposure effect

# Set up a new prediction matrix for the new imputation
ini <- mice(dat_narmice, maxit = 0, print = TRUE)

## Specify the prediction matrix (will include all variables here, as no admin vars and don't want to exclude any)
pred_narmice <- ini$predictorMatrix
pred_narmice

# Make sure the missingness markers don't predict the variables they represent
pred_narmice["infected", "Rtested"] <- 0
pred_narmice

# Set-up predictor matrix for unidentifiable part
# in this case the whole matrix is zeroes because the unidentifiable part of the imputation model contains a single constant (delta*M) rather than additional contributions from the other variables in the dataset
predSens <- ini$predictorMatrix
predSens[predSens == 1] <- 0
predSens

#Set-up list with sensitivity parameter values
pSens <- rep(list(list("")), ncol(dat_narmice))
names(pSens) <- names(dat_narmice)
pSens

# set up vector describing manner of imputation for each variable. As 'infected' is a logistic model and we're adding a sensitivity parameter, we need to specify the 'logregSens' option 
meth_narmice <- ini$method
meth_narmice
meth_narmice["infected"] <- "logregSens"
meth_narmice


# As just 1 variable to impute, can just run one burn-in iteration per imputation
narmice_numimps <- 1
narmice_numiter <- 1

# Number of iterations/simulations to draw bias parameter from
sims <- 10000

# To collect the parameters of interest
res <- as.data.frame(array(dim = c(sims, 22))) # Number of sensitivity values we're going to try, plus the number of parameters we're going to store (here, is 22 [see row below])
colnames(res) <- c("iteration", "csp", "msp", "beta", "se_beta", "sampprev", "beta_sampError", 
                   "cons", "cons_se", "cons_sampError", "men", "men_se", "men_sampError", 
                   "sdage", "sdage_se", "sdage_sampError", "degree", "degree_se", "degree_sampError", 
                   "smoker", "smoker_se", "smoker_sampError")
res

## Set-up prior
prior_mean <- mean_csp
prior_sd <- sd_csp

## Set seed so reproducible
set.seed(842225)

## And time how long this takes
start <- Sys.time()


## Looping over each iteration/simulation, drawing a different bias parameter each time
for (j in 1:sims) {

  # Generate random bias parameter from distribution
  bias <- rnorm(n = 1, mean = prior_mean, sd = prior_sd)
  
  print(paste("Processing iteration", j, "with prior", bias))
  
  # Specify a CSP value for the prediction equation for missing infection data
  pSens[["infected"]]<-list(c(bias))
  
  # NARMICE imputation
  imp_NARMICE <- mice(dat_narmice, m = narmice_numimps, method = meth_narmice, predictorMatrix = pred_narmice,
                      predictorSens=predSens, parmSens=pSens, print = TRUE, maxit = narmice_numiter)
  
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
  mod <- glm(infected.1 ~ sdbmi.1 + men.1 + sdage.1 + degree.1 + smoker.1, 
             family = "binomial", data = imp_NARMICE_wide)
  
  # Store these estimates in the 'results' dataframe
  res[j,"iteration"] <- j
  res[j,"csp"] <- bias
  res[j,"msp"] <- imp_wide_NARMICE_difflogodds
  res[j,"beta"] <- summary(mod)$coefficients["sdbmi.1", 1]
  res[j,"se_beta"] <- summary(mod)$coefficients["sdbmi.1", 2]
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
  res[j,"smoker"] <- summary(mod)$coefficients["smoker.12", 1]
  res[j,"smoker_se"] <- summary(mod)$coefficients["smoker.12", 2]
  res[j,"smoker_sampError"] <- rnorm(n = 1, mean = res[j,"smoker"], sd = res[j,"smoker_se"])
  
  print(res[j,])
  print("")
}

# Check the results table
head(res)

# Calculate time taken
end <- Sys.time()

end - start

## Save these results
write_csv(res, ".\\Results\\NARMICE_Full_results_allIterations.csv")

## Summarise and store these results
res_NARMICE <- as.data.frame(cbind(beta_median = median(res$beta_sampError),
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

res_NARMICE

# Now export results to CSV
write_csv(res_NARMICE, ".\\Results\\results_NARMICE_Full.csv")

    ######################################
    ### Bayesian SM analysis ###

    library(jagsUI)

    set.seed(220511)

    writeLines("
    model
    {
        for (i in 1:N) {
            infected[i] ~ dbern(infected.p[i])
            logit(infected.p[i]) <- beta0 + beta_bmi * sdbmi[i] + 
                beta_men * men[i] + beta_age * sdage[i] + beta_degree * 
                degree[i] + beta_smoker * smoker[i]
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
                psi_bmi * sdbmi[i] + psi_men * men[i] + psi_age * 
                sdage[i] + psi_degree * degree[i] + psi_smoker * 
                smoker[i] + psi_asthma * asthma[i] + psi_diabetes * 
                diabetes[i] + psi_hyperten * hyperten[i]
        }
    # priors
        beta0 ~ dnorm(0.00000E+00, 0.01)
        beta_bmi ~ dnorm(0.00000E+00, 0.01)
        beta_men ~ dnorm(0.00000E+00, 0.01)
        beta_age ~ dnorm(0.00000E+00, 0.01)
        beta_degree ~ dnorm(0.00000E+00, 0.01)
        beta_smoker ~ dnorm(0.00000E+00, 0.01)
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
        psi_men ~ dnorm(0.00000E+00, 0.01)
        psi_age ~ dnorm(0.00000E+00, 0.01)
        psi_degree ~ dnorm(0.00000E+00, 0.01)
        psi_smoker ~ dnorm(0.00000E+00, 0.01)
        psi_asthma ~ dnorm(0.00000E+00, 0.01)
        psi_diabetes ~ dnorm(0.00000E+00, 0.01)
        psi_hyperten ~ dnorm(0.00000E+00, 0.01)
    }
    ", con="/FILEPATH/JAGSmodel.txt")

    ## Set-up prior
    prior_mean <- -mean_csp
    prior_sd <- sd_csp

    dat_SM <- with(dat, list("N"=nrow(dat), "infected"=infected, "sdbmi"=sdbmi, "men"=men, "sdage"=sdage, "degree"=degree, "smoker"=smoker, "asthma"=asthma, "diabetes"=diabetes, "hyperten"=hyperten, "tested"=tested, "prior_mean"=prior_mean, "prior_sd"=prior_sd))
    
    
    # list parmaeters to monitor - just exposure effect and bias parameter?, Dan is saving coefficints in he analysis model and bias parameters only
    Params <- c("beta0", "beta_bmi", "beta_men", "beta_age", "beta_degree", "beta_smoker", "asthma0", "asthma_infected", "asthma_bmi", "asthma_men", "asthma_age", "asthma_degree", "asthma_smoker", "diabetes0", "diabetes_infected", "diabetes_bmi", "diabetes_men", "diabetes_age", "diabetes_degree", "diabetes_smoker", "diabetes_asthma", "hyperten0", "hyperten_infected", "hyperten_bmi", "hyperten_men", "hyperten_age", "hyperten_degree", "hyperten_smoker", "hyperten_asthma", "hyperten_diabetes", "psi0", "psi_infected", "psi_bmi", "psi_men", "psi_age", "psi_degree", "psi_smoker", "psi_asthma", "psi_diabetes", "psi_hyperten")


    # function to generate initial values
    JAGSInits <- function() {list(beta0=rnorm(1, mean=0, sd=1), beta_bmi=rnorm(1, mean=0, sd=1), beta_men=rnorm(1, mean=0, sd=1), beta_age=rnorm(1, mean=0, sd=1), beta_degree=rnorm(1, mean=0, sd=1), beta_smoker=rnorm(1, mean=0, sd=1), asthma0=rnorm(1, mean=0, sd=1), asthma_infected=rnorm(1, mean=0, sd=1), asthma_bmi=rnorm(1, mean=0, sd=1), asthma_men=rnorm(1, mean=0, sd=1), asthma_age=rnorm(1, mean=0, sd=1), asthma_degree=rnorm(1, mean=0, sd=1), asthma_smoker=rnorm(1, mean=0, sd=1), diabetes0=rnorm(1, mean=0, sd=1), diabetes_infected=rnorm(1, mean=0, sd=1), diabetes_bmi=rnorm(1, mean=0, sd=1), diabetes_men=rnorm(1, mean=0, sd=1), diabetes_age=rnorm(1, mean=0, sd=1), diabetes_degree=rnorm(1, mean=0, sd=1), diabetes_smoker=rnorm(1, mean=0, sd=1), diabetes_asthma=rnorm(1, mean=0, sd=1), hyperten0=rnorm(1, mean=0, sd=1), hyperten_infected=rnorm(1, mean=0, sd=1), hyperten_bmi=rnorm(1, mean=0, sd=1), hyperten_men=rnorm(1, mean=0, sd=1), hyperten_age=rnorm(1, mean=0, sd=1), hyperten_degree=rnorm(1, mean=0, sd=1), hyperten_smoker=rnorm(1, mean=0, sd=1), hyperten_asthma=rnorm(1, mean=0, sd=1), hyperten_diabetes=rnorm(1, mean=0, sd=1), psi0=rnorm(1, mean=0, sd=1), psi_infected=rnorm(1, mean=0, sd=1), psi_bmi=rnorm(1, mean=0, sd=1), psi_men=rnorm(1, mean=0, sd=1), psi_age=rnorm(1, mean=0, sd=1), psi_degree=rnorm(1, mean=0, sd=1), psi_smoker=rnorm(1, mean=0, sd=1), psi_asthma=rnorm(1, mean=0, sd=1), psi_diabetes=rnorm(1, mean=0, sd=1), psi_hyperten=rnorm(1, mean=0, sd=1))}

    # Set up number of burnin (5000) and total iteraitons (50000)
    SM_numburnin <- 5000
    SM_numiter <- 50000

    ##### Run JAGS with burn-ins #####
    ## And time how long this takes
    start <- Sys.time()
    res <- jags(data=dat_SM, inits=JAGSInits, parameters.to.save=Params, model.file="/FILEPATH/JAGSModel.txt", n.chains=1, n.iter=SM_numiter, n.burnin=SM_numburnin)

    # Calculate time taken
    end <- Sys.time()

    end - start

    ## Save these results
    write_csv(data.frame(res$samples[[1]]), ".\\Results\\BayesianSM_Full_results_allIterations.csv")

    ## Summarise and store these results
    res_SM <- NULL
    for (param in Params)
    {
    res_SM <- as.data.frame(cbind(res_SM, median=as.vector(unlist(res$q50)[param]), mean=as.vector(unlist(res$mean)[param]), sd=as.vector(unlist(res$sd)[param]), prec_2.5=as.vector(unlist(res$q2.5)[param]), prec_97.5=as.vector(unlist(res$q97.5)[param]), prec_25=as.vector(unlist(res$q25)[param]), prec_75=as.vector(unlist(res$q75)[param])))
    colnames(res_SM)[(-6:0+ncol(res_SM))] <- paste0(param, "_", colnames(res_SM)[(-6:0+ncol(res_SM))])
    }


    # Now export results to CSV
    write_csv(res_SM, ".\\Results\\results_BayesianSM_Full.csv")


