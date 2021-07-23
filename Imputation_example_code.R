# Author: The Dunn Lab, Mass General Hospital, www.thedunnlab.com 
# Description: Script to perform multiple imputation

## imputation for missing exposures and adversities
## Please note: for this example, exposure was a repeated binary measure of maternal psychopathology at 5 separate timepoints and outcome was a continuous measure, SMFQ score at age 10


## load packages ---- 
library(mice) # if this script is causing you trouble, try updating your version of mice
library(tidyr)
library(dplyr)
library(lattice)

## load the imputation function ----

# This function can be set to return 'imp5', which allows you to see if imputation is working correctly based on plots
# When imputation appears successful, 'imp5' should be commented out so that you can create 'imp20' (located below imp5 in the function code), which gives you the final imputed data object
# Please note, while this function performs imputation 20 times, you may want to alter this based on your specific scenario. 
# For a heuristic on how many imputations to perform see: 
# Graham JW, Olchowski AE, Gilreath TD (2007) How many imputations are really needed? Some practical clarifications of multiple imputation theory. Prevention Science 8: 206-213.

source("~/scripts/Imputation/imputation-function-v1.1.R") 

## load the base data for imputation (your original dataset) ---- 
base <- "" # insert R data analytic dataset location here
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
df <- loadRData(base)

## define inputs ----
## confirming that there is no outcome missingness (in example, outcome = SMFQ_10y), if there is, you will subset later to only complete cases
sum(is.na(df$SMFQ_10y))
table(df$SMFQ_10y)

## repeated adversity measure (exposure) - pulls out all columns that begin with the adversity name across the various timepoints. 
# example: mom psychopathology 
mompsy <- grep("^mompsy", names(df), value=T)

# define baseline covariates/confounders
covars <- c("education", "smoking")


## adversity (exposure) of interest
adver <- mompsy
advers <- "mompsy"


## imputing the covariates and exposure variables 

## examine the degree of missingness present in the data frame. Substantial levels of missingness may mean imputation for that variable is not feasible.
missingness <- lapply(covars, function(x){
  N.miss <- sum(is.na(df[,x]))
  perc.miss <- round(N.miss/nrow(df)*100,2)
  var <- x
  vartype <- ifelse(x %in% covars, "covar", "exposure")
  r <- data.frame(cbind(var, N.miss, perc.miss, vartype))
})
#save(missingness, file = "df_covar_missingness_2020.05.13.Rdata")

## imputation variables - ALWAYS include the outcome
# 'imput.vars' are *supplemental* variables used as part of the calculation of the imputed values. 
# They provide extra information for the imputation process. They *are not* imputed for.

imput.vars <- c("SMFQ_10y",
                "alcohol_use",
                "history_of_abuse")

# subset data frame to the variables needed:
df <- df[, colnames(df) %in% c(covars, adver, imput.vars)]
df <- df[complete.cases(df$SMFQ_10y), ] # subset data to those with complete outcome

# You need this adversity list (adv.list) to run the do.impute function
adv.list <- cbind.data.frame(
  adver = advers, # advers previously defined above
 # long.name = advers, # this portion of the command has been deprecated.
  seed1 = c(1234), 
  seed2 = c(5432)
)



#### step 2: perform unstratified imputation ---- 
cont.covars <- c("birthweight") ## list here any continuous variables, either confounders or exposures


## unstratified 
imp.data <- do.impute(df, stratified = FALSE, adver = advers, covars = covars, exp.type = "cat")
#save(imp.data, file = "data/imp20_mompsy_2020.05.18.Rdata")


## option for stratified imputation
## stratified 
imp20.strat <- do.impute(df, stratified = TRUE, adver = advers, covars = covars, exp.type = "cat")
imp20.M <- imp20.strat[[1]] # male
imp20.F <- imp20.strat[[2]] # female
#save(imp20.M, file = "data/imp20_M_20200831.Rdata")
#save(imp20.F, file = "data/imp20_F_20200831.Rdata")



