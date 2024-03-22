### Processing UKBioBank data

## Created 27/1/2022 by Dan Major-Smith
## R version 4.0.4



## Clear workspace, set working directory and load packages
rm(list = ls())

setwd("FILEPATH")

#install.packages("tidyverse")
library(tidyverse)

#install.packages("haven")
library(haven)



## Read in the Stata dataset
dat_raw <- read_dta(file = "FILEPATH\\FILE.dta")

# Make a copy, to save having to read in original again
dat <- dat_raw

summary(dat)


## Prep the variables to match the simulation

# Sex - Rename to 'men', and code women as '0'
table(dat$sex, useNA = "ifany")

dat <- dat %>%
  rename(men = sex) %>%
  mutate(men = ifelse(men == 1, 0, 1))

table(dat$men, useNA = "ifany")

# Age - Standardise by SD - Since start of UKBB (2006)
summary(dat$age_0_0)

summary(dat$age_0_0[!is.na(dat$test_phase1)])

dat$age <- dat$age_0_0

dat <- dat %>%
  rename(sdage = age_0_0) %>%
  mutate(sdage = (sdage - mean(sdage)) / sd(sdage))

summary(dat$sdage)
hist(dat$sdage)

# Degree - recode to 0 = no degree vs 1 = degree
table(dat$eduyears_quali, useNA = "ifany")

dat <- dat %>%
  rename(degree = eduyears_quali) %>%
  mutate(degree = ifelse(degree < 3 & !is.na(degree), 0,
                         ifelse(degree == 3 & !is.na(degree), 1, NA)))

table(dat$degree, useNA = "ifany")

# Smoker - Recode to current smoker or not
table(dat$current_smoke, useNA = "ifany")

dat <- dat %>%
  rename(smoker = current_smoke) %>%
  mutate(smoker = ifelse((smoker == 0 | smoker == 1) & !is.na(smoker), 0,
                         ifelse(smoker == 2 & !is.na(smoker), 1, NA)))

table(dat$smoker, useNA = "ifany")

# BMI - Standardise by SD
summary(dat$bmi_0_0)

summary(dat$bmi_0_0[!is.na(dat$test_phase1)])

dat$bmi <- dat$bmi_0_0

dat <- dat %>%
  rename(sdbmi = bmi_0_0) %>%
  mutate(sdbmi = (sdbmi - mean(sdbmi, na.rm = TRUE)) / sd(sdbmi, na.rm = TRUE))

summary(dat$sdbmi)
hist(dat$sdbmi)

# Infected - Infected with SARS-CoV-2
table(dat$positive_death_negative_phase1, useNA = "ifany")

dat <- dat %>%
  rename(infected = positive_death_negative_phase1)

table(dat$infected, useNA = "ifany")

# Asthma - Diagnosed with asthma
table(dat$resp_diag, useNA = "ifany")

dat <- dat %>%
  rename(asthma = resp_diag)

table(dat$asthma, useNA = "ifany")

# Diabetes - Diagnosed with diabetes
table(dat$diabetes_diag, useNA = "ifany")

dat <- dat %>%
  rename(diabetes = diabetes_diag)

table(dat$diabetes, useNA = "ifany")

# Hypertension - Diagnosed with hypertension
table(dat$hyperten_diag, useNA = "ifany")

dat <- dat %>%
  rename(hyperten = hyperten_diag)

table(dat$hyperten, useNA = "ifany")

# Tested - Tested for SARS-CoV-2 - Are quite a few (25k) missing values. These are all people who died before 2020, so will drop these participants
table(dat$test_phase1, useNA = "ifany")
round((table(dat$test_phase1, useNA = "ifany") / sum(table(dat$test_phase1, useNA = "ifany"))) * 100, 2)

dat <- dat %>%
  rename(tested = test_phase1) %>%
  filter(!is.na(tested))

table(dat$tested, useNA = "ifany")

# Missingness marker for BMI
summary(dat$sdbmi)

dat <- dat %>%
  mutate(Rsdbmi = ifelse(is.na(sdbmi), 1, 0))

table(dat$Rsdbmi, useNA = "ifany")

# Missingness marker for smoker
table(dat$smoker, useNA = "ifany")

dat <- dat %>%
  mutate(Rsmoker = ifelse(is.na(smoker), 1, 0))

table(dat$smoker, useNA = "ifany")

# Missingness marker for degree
table(dat$degree, useNA = "ifany")

dat <- dat %>%
  mutate(Rdegree = ifelse(is.na(degree), 1, 0))

table(dat$Rdegree, useNA = "ifany")

## Drop some unused variables
glimpse(dat)

dat <- dat %>%
  select(-c(n_eid, asthma_SR, diabetes_SR, hcw_testing, positive_death_pop_phase1, death_nonsevere_phase1,
            death_population_phase1, hyperten_base))

summary(dat)

# Check for logical consistency between 'tested' and 'infected' - Are 65 people infected but not tested (these are post-mortem COVID diagnoses. Decided to drop these cases.)
table(dat$infected, dat$tested)
table(dat$infected, dat$tested, useNA = "ifany")

table(dat$infected[dat$tested == 1])
table(dat$infected[dat$tested == 0])

table(dat$tested[dat$infected == 1])
table(dat$tested[dat$infected == 0])

round((table(dat$infected[dat$tested == 0]) / nrow(dat)) * 100, 2)

# Drop these cases
dat <- dat %>%
  mutate(drop = ifelse(infected == 1 & tested == 0, 1, 0)) %>%
  filter(drop == 0 | is.na(drop)) %>%
  select(-drop)

table(dat$infected, dat$tested)
table(dat$tested, useNA = "ifany")
table(dat$infected, useNA = "ifany")


## Summary of variables
table(dat$men, useNA = "ifany")
round(prop.table(table(dat$men)) * 100, 2)

table(dat$smoker, useNA = "ifany")
round(prop.table(table(dat$smoker)) * 100, 2)
round(prop.table(table(dat$smoker, useNA = "ifany")) * 100, 2)

table(dat$degree, useNA = "ifany")
round(prop.table(table(dat$degree)) * 100, 2)
round(prop.table(table(dat$degree, useNA = "ifany")) * 100, 2)

table(dat$asthma, useNA = "ifany")
round(prop.table(table(dat$asthma)) * 100, 2)

table(dat$hyperten, useNA = "ifany")
round(prop.table(table(dat$hyperten)) * 100, 2)

table(dat$diabetes, useNA = "ifany")
round(prop.table(table(dat$diabetes)) * 100, 2)

summary(dat$age)
sd(dat$age)

summary(dat$bmi)
sd(dat$bmi, na.rm = TRUE)
sum(is.na(dat$bmi))
round((sum(is.na(dat$bmi)) / nrow(dat)) * 100, 2)

table(dat$tested, useNA = "ifany")
round(prop.table(table(dat$tested)) * 100, 2)

table(dat$infected, useNA = "ifany")
round(prop.table(table(dat$infected)) * 100, 3)
round(prop.table(table(dat$infected, useNA = "ifany")) * 100, 3)



### Also want to drop anybody with missing covariate data (sdbmi, smoker or degree), and drop associated missingness markers

## Summary of missing data for each covariate

# BMI
sum(is.na(dat$bmi))
round((sum(is.na(dat$bmi)) / nrow(dat)) * 100, 2)

# Smoker
sum(is.na(dat$smoker))
round((sum(is.na(dat$smoker)) / nrow(dat)) * 100, 2)

# Degree
sum(is.na(dat$degree))
round((sum(is.na(dat$degree)) / nrow(dat)) * 100, 2)

# All combined
sum(is.na(dat$bmi) | is.na(dat$smoker) | is.na(dat$degree))
round((sum(is.na(dat$bmi) | is.na(dat$smoker) | is.na(dat$degree)) / nrow(dat)) * 100, 2)

## Remove these cases
dat <- dat %>%
  filter(!is.na(sdbmi) & !is.na(smoker) & !is.na(degree)) %>%
  select(-Rsdbmi, -Rsmoker, -Rdegree)

## Quick summary of data
summary(dat)


## Summary of variables
table(dat$men, useNA = "ifany")
round(prop.table(table(dat$men)) * 100, 2)

table(dat$smoker, useNA = "ifany")
round(prop.table(table(dat$smoker)) * 100, 2)
round(prop.table(table(dat$smoker, useNA = "ifany")) * 100, 2)

table(dat$degree, useNA = "ifany")
round(prop.table(table(dat$degree)) * 100, 2)
round(prop.table(table(dat$degree, useNA = "ifany")) * 100, 2)

table(dat$asthma, useNA = "ifany")
round(prop.table(table(dat$asthma)) * 100, 2)

table(dat$hyperten, useNA = "ifany")
round(prop.table(table(dat$hyperten)) * 100, 2)

table(dat$diabetes, useNA = "ifany")
round(prop.table(table(dat$diabetes)) * 100, 2)

summary(dat$age)
sd(dat$age)

summary(dat$bmi)
sd(dat$bmi, na.rm = TRUE)
sum(is.na(dat$bmi))
round((sum(is.na(dat$bmi)) / nrow(dat)) * 100, 2)

table(dat$tested, useNA = "ifany")
round(prop.table(table(dat$tested)) * 100, 2)

table(dat$infected, useNA = "ifany")
round(prop.table(table(dat$infected)) * 100, 3)
round(prop.table(table(dat$infected, useNA = "ifany")) * 100, 3)


## And now save dataset creating correct z-scores for age and bmi
dat_z <- dat %>%
  select(-c(sdbmi, sdage)) %>%
  rename(sdage = age) %>%
  mutate(sdage = (sdage - mean(sdage)) / sd(sdage)) %>%
  rename(sdbmi = bmi) %>%
  mutate(sdbmi = (sdbmi - mean(sdbmi, na.rm = TRUE)) / sd(sdbmi, na.rm = TRUE))

summary(dat_z)

# Save this data
write_csv(dat_z, "UKB_FullSample_MissingOutcomeOnly_zCorrected.csv")


