
# Author: The Dunn Lab, Mass General Hospital, www.thedunnlab.com 
# Description: Script for a function that runs LARS automatically for imputed data

## select.LARS.imputed: a function that runs LARS automatically for imputed data. This function accomodates sensitive period, accumulation, and recency hypotheses
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


select.LARS.imputed <- function(imp.data, 
                                outcome, 
                                adver, 
                                covars,
                                seed, 
                                hypos, 
                                recency.vec = NULL,
                                inf.method = "fixedLassoInf",
                                exposures = "default",
                                save.derived = FALSE,
                                save.file.name = NULL,
                                step = 1,
                                std.outcome = TRUE) {




  ## Cholesky decomposition functions (averaging the imputed datasets) ----
  bt1 <- expression( {
    # get dataset of relevant variables
    dat <- cbind.data.frame(mget(c(covars, exposures, outcome),inherits=TRUE)) # outcome
    # subset to subject with observed outcome:
    dat <- dat[complete.cases(dat), ]

    # create model matrix (dummy coding factor variables):
    dat2 <- model.matrix(as.formula(paste(outcome, " ~", # outcome
                                          paste(c(covars, exposures),
                                                collapse= "+"))), data=dat)
    # make matrix of variables
    basis <- as.matrix(cbind(dat2[, -1], # remove intercept
                             outcome = dat[,outcome])) # add complete cases outcome

    t(basis) %*% rep(1,n) / n
  } )



  btb <- expression( {
    dat <- cbind.data.frame(mget(c(covars, exposures, outcome), inherits=TRUE)) # outcome 
    dat <- dat[complete.cases(dat), ]

    dat2 <- model.matrix(as.formula(paste(outcome, " ~", # outcome
                                          paste(c(covars, exposures),
                                                collapse= "+"))), data=dat)

    # make matrix of variables
    basis <- as.matrix(cbind(dat2[, -1], # remove intercept
                             outcome = dat[,outcome])) # add complete cases outcome

    t(basis) %*% basis
  } )



  ## data prep ----

  # define sample size:
   n <- sum(!is.na(imp.data$data[,outcome]))

  # define exposures
  if (exposures == "default") {
    exposures <- grep(paste0("^", adver), colnames(imp.data$data), value = T)
    }

  # add ever/never exposed
  # moved the definition of accumulation and recency after decomposition definition
  long <- mice::complete(imp.data, action = "long", include = TRUE)
  exp <- do.call("cbind", lapply(long[, exposures], function(x) as.numeric(as.character(x))))
  long$accumulation <- rowSums(exp, na.rm = F)
  table(long$accumulation)
  # long %>% filter(is.na(accumulation)) %>% group_by(.imp) %>% summarise(N=n())
  
  # confirming that all missingness occurred in the original (unimputed) data frame
  exp.sum1 <- rowSums(exp, na.rm = T)
  long$ever <- ifelse(exp.sum1 > 0, 1, ifelse(long$accumulation == 0, 0, NA))

  if ("recency" %in% hypos){
    long$recency <-  rowSums(exp %*% diag(recency.vec))
  } else {
    # if you do not need to test for recency, force the entire vector to be a normal random variable
    long$recency <- rnorm(n = nrow(long))
  }

  if (std.outcome == TRUE) {
    long[,outcome] <- (long[,outcome] - mean(long[,outcome]))/sd(long[,outcome])
  }

  # summary(long$recency)
  imp.data <- as.mids(long)  # you must run this line after creating new variables in the imp.data long dataset

  if (save.derived == TRUE){
    save(imp.data, file = save.file.name)
  }

  # Cholesky decomposition
  bt1s <- with(data=imp.data, eval(bt1))
  btbs <- with(data=imp.data, eval(btb))
  Vmat <- apply(simplify2array(btbs$analyses),1:2,mean)
  Mvec <- apply(simplify2array(bt1s$analyses),1:2,mean) 
  set.seed(seed)
  z <- rnorm(n)
  Z <- lm(z ~ 1)$residuals
  Z <- Z / sqrt(sum(Z^2))
  for(j in 2:length(Mvec)) {
    z <- rnorm(n)
    Z <- cbind(Z, lm(z ~ Z)$residuals)
    Z[,j] <- Z[,j] / sqrt(sum(Z[,j]^2))
  }

  Amat <- chol(Vmat - n * Mvec %*% t(Mvec))
  New <- Z %*% Amat + rep(1,n) %*% t(Mvec)
  New <- as.data.frame(New)
  

 exp <- do.call("cbind", lapply(New[, paste0(exposures, "1")], function(x) as.numeric(as.character(x))))
 New$accumulation <-  rowSums(exp, na.rm = F)
  if ("recency" %in% hypos){
    New$recency <-  rowSums(exp %*% diag(recency.vec), na.rm = F)
  } else {
    # if you do not need to test for recency, force the entire vector to be a normal random variable
   New$recency <- rnorm(n = nrow(New))
  }


  ## select covariates -- note that the data frame New should only have covariates, exposures, and outcome
  covariates.names <- setdiff(colnames(New), c(paste0(exposures, "1"),
                                               "outcome", "accumulation", "recency"))
  covariates <- as.matrix(New[,covariates.names])

  X_hypos <- as.matrix(New[, c(paste0(exposures, "1"), hypos)])
  y <- New$outcome
  X_residual <- lm(X_hypos ~ covariates)$resid
  
  ## Frish-Waugh-Lovell - before lasso
  y <- lm(y~covariates)$resid

  lasso <- lars(X_residual, y)



  ## model selection ----
  
  # creation of elbow plots
  
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
  plot(c(0, variables), lasso$R2, type='l',
       xlab="Variables selected", ylab="R-squared")
  text(c(0, variables), rep(0, max(variables)+1),
       labels=c("", additions), srt=90, adj=0)
  

  ## fixed Lasso inference ----
  if ("fixedLassoInf" %in% inf.method) {


    sumsq <- lasso$normx
    X_normed <- scale(X_residual, scale = sumsq)

    ## include more steps
    fli <- fixedLassoInf(X_normed, y, lasso$beta[step+1,], lasso$lambda[step+1], alpha = 0.05)
    scale <- sumsq[fli$vars]
    fli$coef0 <- fli$coef0/scale
    fli$sd <- fli$sd/scale
    fli$ci <- fli$ci/cbind(scale,scale)
    fli$vars <- names(fli$vars)
    print("selectiveInference result:")
    fli.res <- as.data.frame(cbind(fli$vars, fli$coef0, fli$sd, fli$pv, fli$ci))
    colnames(fli.res) <- c("Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")

    fli.res[,-1] <- lapply(fli.res[,-1], function(x) round(as.numeric(as.character(x)), 4))
    print(fli.res)
    print(knitr::kable(fli.res))

    print("p-value with more decimal points:")
    print(fli$pv)

    
  }

}
