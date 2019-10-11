
## example showing the use of the LARS function
## for complete case analysis 



# simulate toy data -------------------------------------------------------


n <- 2000

set.seed(42)


# simulate exposure data 
abuse_18m <- rbinom(n, 1, 0.17)
odds <- 0.1115/(1-0.1115) * 8^abuse_18m
abuse_30m <- rbinom(n, 1, odds/(1+odds))
odds <- 0.111/(1-0.111) * 8^abuse_30m 
abuse_42m <- rbinom(n, 1, odds/(1+odds))
odds <- 0.111/(1-0.111) * 8^abuse_42m
abuse_57m <- rbinom(n, 1, odds/(1+odds))
odds <- 0.111/(1-0.111) * 8^abuse_57m
abuse_69m <- rbinom(n, 1, odds/(1+odds))  


# simulate covariates 
education <- rbinom(n, 3, 0.2)
smoking <- rbinom(n, 2, 0.09)

# simulate outcome data 
# associated with the first exposure 
y <- abuse_18m*0.2 + education*0.05 + smoking*0.1 + rnorm(n) 

dat <- as.data.frame(cbind(abuse_18m, abuse_30m, abuse_42m, abuse_57m, abuse_69m, 
             education, smoking, y))
dat$education <- as.factor(dat$education)
dat$smoking <- as.factor(dat$smoking)


# source LARS function ----------------------------------------------------
source("function-SLCMA-complete-cases.R")



# run analysis ------------------------------------------------------------

# note that the code was designed to accomodate high dimensional outcomes 
# therefore the format is a bit odd if you only have a single outcome

results <- select.LARS.complete(dat = dat, # dataset 
                     outcome.vec = c("y"), # name(s) of your outcome measure
                     adver = "abuse",  # format names of your exposures so that they start with this string 
                     covars = c("education", "smoking"), # covariates 
                     hypos = c("accumulation", "recency"), # additive hypotheses  
                     recency.vec = c(18/12, 30/12, 42/12, 57/12, 69/12), # weights for recency 
                     inf.method = "covTest", # inference method, covTest is default 
                    exposures = "default") # leave as default 


