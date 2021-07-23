
# Author: The Dunn Lab, Mass General Hospital, www.thedunnlab.com 
# Description: Script to execute SLCMA function on complete case data. This script focuses on a scenario with one category of exposure.
# No source function for this script as function included within the script

## Set up
# Load libraries
library(lars)
library(tidyr)
library(knitr)
library(dplyr)
library(selectiveInference)

# load data
# a function that lets you rename the object loaded
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

df <- loadRData("") # here, provide file path to your data frame


# defining criteria for fn
outcome = "" # character string of name for outcome variable, e.g., "y" or "depression_score", must be same name as the respective column in your data frame.
adver = ""  # name for exposure points of interest, for example, if you have 5 timepoints where exposure to abuse was measured 
            # they should have the same naming convention, e.g., "abuse_18mo", "abuse_30mo", 
            # then write "abuse" here (with quotes).

exposures <- grep(paste0("^", adver), colnames(df), value = T) # this command creates a vector listing all sensitive period measurements containing "abuse" in the name (if following example above)
covars = c("education", "smoking") # list of covariates of interest
hypos = c("accumulation", "recency") # life course hypotheses you'd like to consider besides sensitive periods
recency.vec = list(c(18, 30, 42, 57, 69)/12) # list of timepoints when the exposure was measured to create the recency hypothesis. 
                                            # As an example here, exposure was measured at 18 months, 30 months, etc. and will be evaluated in years, 
                                            # thus, all values are divided by 12. Create a vector for each adversity studied.
names(recency.vec) <- adver
step = 1 # number of hypotheses you will accept (this can be changed based on the results of the elbow plot)

#define df
df <- df %>% filter(rowSums(is.na(.[,c(exposures,covars)])) == 0) # filtering dataset to the complete cases

# add accumulation and recency hypotheses to df
exp <- do.call("cbind", lapply(df[, exposures], function(x) as.numeric(as.character(x))))
df <- df %>% 
  mutate(accumulation = rowSums(exp)) %>% 
  mutate(recency = rowSums(exp %*% diag(recency.vec)))


#### NOTE: Typical user should not have to change any code beyond this point in the script all that is required is to run all remaining lines now that the inputs above have been assigned --------
## select covariates -- 
## what are the categorical covariates - your variables should be categorized appropriately prior to SLCMA analysis
cat.covars <- do.call("c", lapply(1:length(covars), function(x) {
  cov1 <- covars[x]
  if (is.factor(df[, cov1])) {
    return(cov1)
  } else {
    return(NULL)
  }
}))


cat.covariates <- model.matrix(as.formula(paste0("~", paste(cat.covars, collapse = "+"))), df)[,-1]
covariates <- as.matrix(cbind(cat.covariates, as.matrix(df[,setdiff(covars, cat.covars)])))
X_hypos <- do.call("cbind", lapply(df[, c(exposures, hypos)], function(x) as.numeric(as.character(x))))
X_residual <- lm(X_hypos ~ covariates)$resid

y <- df[, outcome]
y <- as.numeric(as.character(y))

## Frish-Waugh-Lovell - apply before lasso - reduces confounding bias
y <- lm(y~covariates)$resid

lasso <- lars(X_residual, y) # penalization method
last.action <- length(lasso$actions)
variables <- numeric(last.action)
selection <- vector("list", last.action)

for(action in 1:last.action) {
  variables[action] <- sum(lasso$beta[as.character(action), ] != 0)
  selection[[action]] <- attributes(which(lasso$beta[as.character(action),
                                                     ] != 0))$names
}

additions <- character(last.action)
for(action in 2:last.action) {
  current_selection <- lasso$beta[as.character(action-1), ] != 0
  new_selection <- lasso$beta[as.character(action), ] != 0
  if(variables[action] > variables[action-1]) {
    additions[action] <-
      dimnames(X_residual)[[2]][new_selection != current_selection]
  }
}

additions[1] <- selection[1]

# Elbow plot
par(mar=c(5,4,4,5)+0.1)
plot(c(0, variables), lasso$R2, type='l', # lasso$R2 - elbow plots model R2 value
     xlab="Variables selected", ylab="R-squared")
text(c(0, variables), rep(0, max(variables)+1),
     labels=c("", additions), srt=90, adj=0)


# sI - post-selection inference method
sumsq <- lasso$normx
X_normed <- scale(X_residual, scale = sumsq)

# creating output
fli <- fixedLassoInf(X_normed, y, lasso$beta[step+1,], lasso$lambda[step+1], alpha = 0.05)
scale <- sumsq[fli$vars]
fli$coef0 <- fli$coef0/scale
fli$sd <- fli$sd/scale
fli$ci <- fli$ci/cbind(scale,scale)
fli$vars <- names(fli$vars)
print("selectiveInference result:")
#print(fli)
fli.res <- as.data.frame(cbind(fli$vars, fli$coef0, fli$sd, fli$pv, fli$ci))
colnames(fli.res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

fli.res[,-1] <- lapply(fli.res[,-1], function(x) round(as.numeric(as.character(x)), 4))

print(fli.res)
sI.res <- knitr::kable(fli.res)
print(knitr::kable(fli.res))

print("p-value with more decimal points:")
print(fli$pv)    
