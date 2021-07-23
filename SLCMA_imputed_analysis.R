
# Author: The Dunn Lab, Mass General Hospital, www.thedunnlab.com 
# Description: Script to execute SLCMA function on imputed data. This script focuses on a scenario with one category of exposure.

# Description of function inputs:
## imp.data: multiply imputed dataset using MICE
## outcome: outcome variable name (string)
## adver: short name of adversity (string)
## covars: a list of covariates
## seed: random seed
## hypos: additional hypotheses to be tested (aside from the individual time points), e.g., accumulation hypothesis
## recency.vec: recency vector - weighting by age (in years, for our example)
## inf.method: post-selective inference method to be used after LARS. We use selective Inference, denoted here as "fixedLassoInf"
## exposures: default to all individual time points; otherwise specify a list of collapsed time points
## save.derived: whether to save the derived imp.data object for subsequent analyses
## save.file.name: if saving derived imp.data, what's the name of the file
## step: how many steps you want to estimate, in other words, how many key variables do you want to include in final model? This can be altered based on the results of the elbow plot.
## std.outcome: whether the outcome should be standardized



## Set up
# Load libraries
library(mice)
library(lars)
library(tidyr)
library(knitr)
library(dplyr)
library(selectiveInference)


## Load data source

# a function that lets you rename the object loaded
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

imp.data <- loadRData("") # here, provide file path to your MICE imputed data object


## Source LARS function
source("~/LARS-function-imputed_data.R") # provide file path to LARS function script available on Dunn Lab's GitHub page

# Load inputs
outcome = "" # character string of name for outcome variable, e.g., "y" or "depression_score", must be same name as in your data frame.
adver = ""  # name for exposure points of interest, for example, if you have 5 timepoints where exposure to abuse was measured 
            # and they have the same naming convention, e.g., "abuse_18mo", "abuse_30mo", 
            # you would store these names in a character vector named abuse (no quotes) and write "abuse" here (with quotes).

# defining recency hyp parameters, one vector for each adversity being studied
recency.vec <- list(c(18, 30, 42, 57, 69)/12) # provide list of timepoints when the exposure was measured to create the recency hypothesis. 
                                              # As an example here, exposure was measured at 18 months, 30 months, etc. and will be evaluated in years, 
                                              # thus, all values are divided by 12.
names(recency.vec) <- adver


## imputed N: 
results_lars <- select.LARS.imputed(imp.data = imp.data, 
                                    outcome = outcome, 
                                    adver = adver, 
                                    exposures = "default",
                                    recency.vec = recency.vec,
                                    covars = c("education", "smoking"), # list covariates here
                                    seed = 1234, 
                                    hypos = c("accumulation", "recency"), 
                                    inf.method = "fixedLassoInf",
                                    save.derived = FALSE,
                                    save.file.name = NULL,
                                    step = 1,
                                    std.outcome = FALSE)


