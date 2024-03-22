### This script shows how to apply the Bayesian SM analysis to the simulated data as described in our paper.

## Created by Emily Kawabata
## R version 4.1.0

## Clear workspace, set working directory and load packages
rm(list=ls())

# Set working directory
setwd("FILEPATH")

# Load in R packages
#install.packages("tidyverse") 
library(tidyverse)
#install.packages("jagsUI")
library(jagsUI)

# for the ith simulated data set: i=1,..., 1000
i <- 1

# Read in the full dataset with no missing covariates - this will differ, depending on whether the selection model or pattern-mixture model simulations are used, and whether the causal or null datasets are used
data <- read_csv(paste0("100k_causal_csv/Dataset_", i, ".csv"))
#data <- read_csv(paste0("100k_null_csv/Dataset_", i, ".csv"))
#data <- read_csv(paste0("PMM_100k_causal_csv/Dataset_", i, ".csv"))
#data <- read_csv(paste0("PMM_100k_null_csv/Dataset_", i, ".csv"))

# Make a variable for infection status which is missing if not tested
data$infected[data$itest == 0] <- NA

## Set the seed - Using random seeds generated in Stata for this (each simulation will have a different set of seeds)
df_seed <- read_csv("SEED_FILE_SMCausalVeryInf.csv")
#df_seed <- read_csv("SEED_FILE_SMCausalInf.csv")
#df_seed <- read_csv("SEED_FILE_SMCausalVague.csv")
#df_seed <- read_csv("SEED_FILE_SMNullVeryInf.csv")
#df_seed <- read_csv("SEED_FILE_SMNullInf.csv")
#df_seed <- read_csv("SEED_FILE_SMNullVague.csv")
#df_seed <- read_csv("SEED_FILE_PMMCausalVeryInf.csv")
#df_seed <- read_csv("SEED_FILE_PMMCausalInf.csv")
#df_seed <- read_csv("SEED_FILE_PMMCausalVague.csv")
#df_seed <- read_csv("SEED_FILE_PMMNullVeryInf.csv")
#df_seed <- read_csv("SEED_FILE_PMMNullInf.csv")
#df_seed <- read_csv("SEED_FILE_PMMNullVeryInf.csv")

# Set seed number
seed <- as.numeric(df_seed[i, ])
print(seed)
set.seed(seed)
rm(df_seed)

## Set up to run JAGS
## Set up hyperparamters for the bias parameter - Select either very informative prior (mean = 7.849, SD = 1), informative prior (mean = 7.849, SD = 2) or vague prior (mean = 0, SD = 10)
prior_mean <- 7.849
#prior_mean <- 0

prior_sd <- 1
#prior_sd <- 2
#prior_sd <- 10

# List data to be used in JAGS
data_SM <- with(data, list("N"=nrow(data), "infected"=infected, "sdbmi"=incomplete_sdbmi, "men"=men, "sdage"=sdage, "degree"=degree, "smoker"=incomplete_smoker, "asthma"=asthma, "diabetes"=diabetes, "hyperten"=incomplete_hyperten, "tested"=itest, "prior_mean"=prior_mean, "prior_sd"=prior_sd))

# List functions to generate initial values
JAGSInits <- function() {list(beta0=rnorm(1, mean=0, sd=1), beta_bmi=rnorm(1, mean=0, sd=1), beta_men=rnorm(1, mean=0, sd=1), beta_age=rnorm(1, mean=0, sd=1), beta_degree=rnorm(1, mean=0, sd=1), beta_smoker=rnorm(1, mean=0, sd=1), bmi0=rnorm(1, mean=0, sd=1), bmi_men=rnorm(1, mean=0, sd=1), bmi_age=rnorm(1, mean=0, sd=1), bmi_degree=rnorm(1, mean=0, sd=1), bmi_smoker=rnorm(1, mean=0, sd=1), bmi_prec=rgamma(1, shape=0.01, rate=0.01), smoker0=rnorm(1, mean=0, sd=1), smoker_men=rnorm(1, mean=0, sd=1), smoker_age=rnorm(1, mean=0, sd=1), smoker_degree=rnorm(1, mean=0, sd=1), asthma0=rnorm(1, mean=0, sd=1), asthma_infected=rnorm(1, mean=0, sd=1), asthma_bmi=rnorm(1, mean=0, sd=1), asthma_men=rnorm(1, mean=0, sd=1), asthma_age=rnorm(1, mean=0, sd=1), asthma_degree=rnorm(1, mean=0, sd=1), asthma_smoker=rnorm(1, mean=0, sd=1), diabetes0=rnorm(1, mean=0, sd=1), diabetes_infected=rnorm(1, mean=0, sd=1), diabetes_bmi=rnorm(1, mean=0, sd=1), diabetes_men=rnorm(1, mean=0, sd=1), diabetes_age=rnorm(1, mean=0, sd=1), diabetes_degree=rnorm(1, mean=0, sd=1), diabetes_smoker=rnorm(1, mean=0, sd=1), diabetes_asthma=rnorm(1, mean=0, sd=1), hyperten0=rnorm(1, mean=0, sd=1), hyperten_infected=rnorm(1, mean=0, sd=1), hyperten_bmi=rnorm(1, mean=0, sd=1), hyperten_men=rnorm(1, mean=0, sd=1), hyperten_age=rnorm(1, mean=0, sd=1), hyperten_degree=rnorm(1, mean=0, sd=1), hyperten_smoker=rnorm(1, mean=0, sd=1), hyperten_asthma=rnorm(1, mean=0, sd=1), hyperten_diabetes=rnorm(1, mean=0, sd=1), psi0=rnorm(1, mean=0, sd=1), psi_infected=rnorm(1, mean=0, sd=1), psi_bmi=rnorm(1, mean=0, sd=1), psi_asthma=rnorm(1, mean=0, sd=1), psi_diabetes=rnorm(1, mean=0, sd=1), psi_hyperten=rnorm(1, mean=0, sd=1))}

# List parameters to monitor
Params <- c("beta0", "beta_bmi", "beta_men", "beta_age", "beta_degree", "beta_smoker", "asthma0", "asthma_infected", "asthma_bmi", "asthma_men", "asthma_age", "asthma_degree", "asthma_smoker", "diabetes0", "diabetes_infected", "diabetes_bmi", "diabetes_men", "diabetes_age", "diabetes_degree", "diabetes_smoker", "diabetes_asthma", "hyperten0", "hyperten_infected", "hyperten_bmi", "hyperten_men", "hyperten_age", "hyperten_degree", "hyperten_smoker", "hyperten_asthma", "hyperten_diabetes", "psi0", "psi_infected", "psi_bmi", "psi_asthma", "psi_diabetes", "psi_hyperten")

# Create a model to be fitted in JAGS
writeLines("
model
{
    for (i in 1:N) {
        infected[i] ~ dbern(infected.p[i])
        logit(infected.p[i]) <- beta0 + beta_bmi * sdbmi[i] + beta_men * men[i] + beta_age * sdage[i] + beta_degree * degree[i] + beta_smoker * smoker[i]
        sdbmi[i] ~ dnorm(bmi0 + bmi_men * men[i] + bmi_age * sdage[i] + bmi_degree * degree[i] + bmi_smoker * smoker[i], bmi_prec)
        smoker[i] ~ dbern(smoker.p[i])
        logit(smoker.p[i]) <- smoker0 + smoker_men * men[i] + smoker_age * sdage[i] + smoker_degree * degree[i]
        asthma[i] ~ dbern(asthma.p[i])
        logit(asthma.p[i]) <- asthma0 + asthma_infected * infected[i] + asthma_bmi * sdbmi[i] + asthma_men * men[i] + asthma_age * sdage[i] + asthma_degree * degree[i] + asthma_smoker * smoker[i]
        diabetes[i] ~ dbern(diabetes.p[i])
        logit(diabetes.p[i]) <- diabetes0 + diabetes_infected * infected[i] + diabetes_bmi * sdbmi[i] + diabetes_men * men[i] + diabetes_age * sdage[i] + diabetes_degree * degree[i] + diabetes_smoker * smoker[i] + diabetes_asthma * asthma[i]
        hyperten[i] ~ dbern(hyperten.p[i])
        logit(hyperten.p[i]) <- hyperten0 + hyperten_infected * infected[i] + hyperten_bmi * sdbmi[i] + hyperten_men * men[i] + hyperten_age * sdage[i] + hyperten_degree * degree[i] + hyperten_smoker * smoker[i] + hyperten_asthma * asthma[i] + hyperten_diabetes * diabetes[i]
        tested[i] ~ dbern(tested.p[i])
        logit(tested.p[i]) <- psi0 + psi_infected * infected[i] + psi_bmi * sdbmi[i] + psi_asthma * asthma[i] + psi_diabetes * diabetes[i] + psi_hyperten * hyperten[i]
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
", con="SimulationJAGSmodel.txt")

# Set up number of burn-in (5000) and total iterations (50000)
SM_numburnin <- 5000
SM_numiter <- 50000

## Run JAGS
start <- Sys.time()

# Run JAGS - this may take hours to run
res <- jags(data=data_SM, inits=JAGSInits, parameters.to.save=Params, model.file="SimulationJAGSmodel.txt", n.chains=1, n.iter=SM_numiter, n.burnin=SM_numburnin, verbose=TRUE)

# Calculate time taken
end <- Sys.time()
end - start

# Save results
## Summarise results
res_SM <- NULL
for (param in Params)
{
    res_SM <- as.data.frame(cbind(res_SM, median=as.vector(unlist(res$q50)[param]), mean=as.vector(unlist(res$mean)[param]), sd=as.vector(unlist(res$sd)[param]), prec_2.5=as.vector(unlist(res$q2.5)[param]), prec_97.5=as.vector(unlist(res$q97.5)[param]), prec_25=as.vector(unlist(res$q25)[param]), prec_75=as.vector(unlist(res$q75)[param])))
    colnames(res_SM)[(-6:0+ncol(res_SM))] <- paste0(param, "_", colnames(res_SM)[(-6:0+ncol(res_SM))])
}

# Save results
write_csv(res_SM, paste0("SMresults_", i, ".csv"))