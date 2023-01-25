
# Author: The Dunn Lab, Mass General Hospital, www.thedunnlab.com 
# Description: Script for a function that has stratified and unstratified imputation methods included

## df: whole analytic dataset 
## adver: adversity (repeated exposure) that's being imputed on
## covars: a list of covariates 
## exp.type: type of exposure, default = categorical 



# v1.1: added plotting to visually examine the imputation


do.impute <- function(df, stratified, adver, covars, exp.type, env = globalenv()){
  
  ## impute separately if stratified by sex (or a different covariate of interest)
  if (stratified == "TRUE"){
    
    covars <- covars[covars != "Female"]
    
    df1 <- df %>% filter(Female == 0)
    imp20.male <- do.call(impute, list(df1, adver, covars, exp.type))
    
    df2 <- df %>% filter(Female == 1)
    imp20.female <- do.call(impute, list(df2, adver, covars, exp.type))
    
    
    return(list(imp20.male, imp20.female))
  } else {
    imp20.unstrat <- do.call(impute, list(df, adver, covars, exp.type))
    return(imp20.unstrat)
  }
}


## a function called in the do.impute function ----
## performs imputation on one data frame 

impute <- function(df, adver, covars, exp.type, env = globalenv()){
  
  # turn categorical variables into factors - define cat vars then convert to factor
  # better safe than sorry, recode your variables first in your analytic dataset! 
  if (exp.type == "cat"){
    cat.vars <- c(setdiff(covars, cont.covars), unlist(mget(advers, envir = globalenv())))
  } else {
    cat.vars <- c(setdiff(covars, cont.covars))
  }
  df[,cat.vars] <- lapply(df[,cat.vars], as.factor)
 
  # set method for each exposure with missing data:
  ini <- mice(df, max = 0, print = FALSE)
  meth <- ini$meth
  meth[] <- "" 
  exposures <- get(adver, envir = globalenv()) #exposures becomes a vector string
  # which exposures having missing data?
  print(exposures.miss <- names(which(apply(df[, exposures], 2, 
                                            function(j) sum(is.na(j))) > 0))) # pulls columns in df using these strings
  
  ## ask the user whether the exposure variable is categorical 

  if (exp.type == "cat"){
    ## specify the imputation method for exposure variables 
    var.length <- do.call(rbind, lapply(df[, exposures.miss],function(j) length(table(j))))
    meth[exposures.miss[var.length == 2]] <- "logreg"
    meth[exposures.miss[var.length > 2]] <- "polyreg" 
  } else {
    meth[exposures.miss] <- "pmm"
  }

  print(covars.miss <- names(which(apply(df[, covars], 2, 
                                         function(j) sum(is.na(j))) > 0)))
  # how many levels do the covariates with missing data have?
  covars.miss.lngth <- do.call(rbind, lapply(apply(df[, covars.miss], 2, table), 
                                             function(j) length(j)))
  
  # set imputation method for covars with missing data with two levels to 
  # logreg, and for covars with >2 levels to polyreg, and any continuous 
  # covariates to pmm:
  meth[covars.miss[covars.miss.lngth == 2]] <- "logreg"
  meth[covars.miss[covars.miss.lngth > 2]] <- "polyreg"
  
  for (cont.covar in cont.covars){
    meth[cont.covar] <- "pmm"
  }
  meth
  # create predictor matrix at default minimum bivariate correlation (r > 0.1)
  pred <- quickpred(df)

# run initial imputation of 5 datasets for 50 iterations and check 
# distribution of imputed data.
 imp5 <- mice(df, meth = meth, pred = pred, m = 5, maxit = 50, 
              seed = adv.list[adv.list$adver == adver, "seed1"], print=F)
   
  # quick visual examination of the imputed values
# print(plot(imp5, exposures.miss, 
#       main="Imputation at min. corr. = 0.1 (default)"))
# print(plot(imp5, covars.miss, 
#         main="Imputation at min. corr. = 0.1 (default)"))
  
# print(densityplot(x=imp5, data=as.formula(paste("~", paste(exposures.miss, collapse= "+")))))
# print(densityplot(x=imp5, data=as.formula(paste("~", paste(covars.miss, collapse= "+")))))
#   
# rm(imp5, ini)
  
  # once satisfied with the outcome of imp5, run true imputation of 20 datasets at 25 iterations and save each dataset:
#imp20 <- mice(df, meth = meth, pred = pred, m = 20, maxit = 25, 
  #   seed = adv.list[adv.list$adver == adver, "seed2"], print=F)
#return(imp20)
}
