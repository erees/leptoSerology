library(Rsero)
library(Rcpp)
library(tidyverse)
source("code/models/timeVaryingFOIFunction.R")

############################################################################
# Read in and process full data
# Consists of line level data with serostatus, age, and division
## Data not available for sharing
############################################################################

# Get the data in the appropriate format for RSero package
age <- seroDat$age
seropositive <- seroDat$seropositive
seroModelData <- SeroData(age_at_sampling = age, Y = seropositive, sampling_year = 2013) 

seroprevalence(seroModelData)

#####################################################################
# Comparing catalytic and reverse catalytic models to Rjags code
#####################################################################

#####################################################################
## 1. Catalytic model
#####################################################################

## Define model (Using default priors for FOI)
ConstantModel = FOImodel(type = 'constant') 

FOIfit.constant = fit(data = seroModelData, ## specify data 
                      model = ConstantModel, ## specify model
                      chains=4, ## number of chains
                      iter = 10000) ## numer of iterations

seroprevalence.fit(YLIM=0.8, FOIfit = c(FOIfit.constant, FOIfit.constant)) ## plot output
parameters_credible_intervals(FOIfit.constant) ## print paramaters

## Using this to check for convergence of chains, neff Rhat, and posterior density
stanObject <- FOIfit.constant$fit
launch_shinystan(stanObject)
stan_dens(stanObject, pars = c("foi[1]"))

logLikelihood <- extractLogLik(FOIfit.constant)

loo::waic(logLikelihood)
loo(logLikelihood)

## Saving the model results
saveRDS(FOIfit.constant, "stan_constant")
FOIfit.constant <- readRDS("stan_constant")

#####################################################################
## 2. Reverse catalytic model
#####################################################################

## Define model (Using default priors for FOI and Rho [waning])
ConstantModelReverse = FOImodel(type = 'constant', ## type of model
                                seroreversion = 1) ## seroreversion is TRUE

FOIfit.constantReverse = fit(data = seroModelData,  
                             model = ConstantModelReverse, 
                             chains=4,
                             iter = 10000)

seroprevalence.fit(FOIfit.constantReverse, YLIM=0.6)
parameters_credible_intervals(FOIfit.constantReverse)

## Using this to check for convergence of chains, neff Rhat, and posterior density
stanObject <- FOIfit.constantReverse$fit
launch_shinystan(stanObject)
pairs(FOIfit.constantReverse$fit, pars = c("foi[1]", "rho"))
stan_dens(stanObject, pars = c("foi[1]","rho"))

## Saving the model results
saveRDS(FOIfit.constantReverse, "savedModels/stan_constant_reverse")
FOIfit.constantReverse <- readRDS("savedModels/stan_constant_reverse")

## Extract log liklihood to calculate WAIC
logLikelihood <- extractLogLik(FOIfit.constantReverse)

loo::waic(logLikelihood)


# Comapre the catalytic and reverse catalytic using AIC and DIC
DIC.Reverse = compute_information_criteria(FOIfit.constantReverse)
DIC.constant = compute_information_criteria(FOIfit.constant)



#####################################################################
# Allowing for time-varying FOI (allowing for outbreaks to occur)
# Comparing three models which include an outbreak
# 3. Constant FOI with 1 outbreak in last 2 years
# 4. Constant FOI with 1 outbreak in last 5 years
# 5. No constant FOI but with 1 outbreak in last 10 years
#####################################################################


#####################################################################
# 3. Constant FOI with 1 outbreak in last 2 years
#####################################################################

## Define model and specify priors
contantOutbreakModel1 = FOImodel(type='constantoutbreak', ## Type of model 
                                 K=1, ## number of outbreaks
                                 seroreversion = 1, ## seroconversion is TRUE
                                 priorT1 = 0, ## lower bound for outbreak time (time since sampling survey) (uniform distribution)
                                 priorT2 = 2, ## upper bound for outbreak time (time since sampling survey) (uniform distribution)
                                 priorRho1 = 0, ## lower bound for waning estimate (uniform distribution)
                                 priorRho2 = 10) ## upper bound for waning estimate (uniform distribution)

FOIfit.constantOutbreak1 = fit(data = seroModelData, ## specify model
                               model = contantOutbreakModel1, ## model type
                               chains = 4, # number of chains
                               iter = 100000, # number of itterations
                               thin = 10, 
                               warmup = 20000
)


parameters_credible_intervals(FOIfit.constantOutbreak1) # get parameters
# Calculate AIC, DIC and WAIC
DIC.constantOutbreak1 = compute_information_criteria(FOIfit.constantOutbreak1)

seroprevalence.fit(FOIfit.constantOutbreak1, YLIM = 0.8) # plot model

FOIfit.constantOutbreak1$fit
pairs(FOIfit.constantOutbreak1$fit, pars = c("T[1]", "constant", "rho"))


stanObject_constantOutbreak1 <- FOIfit.constantOutbreak1$fit
fit_summary <- summary(stanObject_constantOutbreak1)
test <- fit_summary$summary
fit_summary$c_summary

## Check for convergence of chains, neff, Rhat and posteriors of all parameters
# Want ESS to be >200 
# Want Rhat <1.1
stan_diag(stanObject_constantOutbreak1)
stan_trace(stanObject_constantOutbreak1,pars = "T[1]")
stan_trace(stanObject_constantOutbreak1,pars = "rho")
stan_trace(stanObject_constantOutbreak1,pars = "constant")
stan_dens(stanObject_constantOutbreak1, pars = c("T[1]","rho","constant"))

launch_shinystan(stanObject_constantOutbreak1)

# save model
saveRDS(FOIfit.constantOutbreak1, "stan_constantOutbreak1_reverse_2_years")

FOIfit.constantOutbreak1 <- readRDS("savedModels/stan_constantOutbreak1_reverse_2_years")

## Extract log liklihood to calculate WAIC
logLikelihood_constantOutbreak1 <- extractLogLik(FOIfit.constantOutbreak1)

loo::waic(logLikelihood_constantOutbreak1)

#####################################################################
# 4. Constant FOI with 1 outbreak in last 5 years
#####################################################################

contantOutbreakModel2 = FOImodel(type='constantoutbreak', ## Type of model 
                                 K=1, ## number of outbreaks
                                 seroreversion = 1, ## seroconversion is TRUE
                                 priorT1 = 0, ## lower bound for outbreak time (time since sampling survey) (uniform distribution)
                                 priorT2 = 5, ## upper bound for outbreak time (time since sampling survey) (uniform distribution)
                                 priorRho1 = 0, ## lower bound for waning estimate (uniform distribution)
                                 priorRho2 = 10) ## upper bound for waning estimate (uniform distribution)

FOIfit.constantOutbreak2 = fit(data = seroModelData, ## specify model
                               model = contantOutbreakModel2, ## model type
                               chains = 4, # number of chains
                               iter = 100000, # number of itterations
                               thin = 10, 
                               warmup = 20000
)


parameters_credible_intervals(FOIfit.constantOutbreak2) # get parameters
# Calculate AIC, DIC and WAIC
DIC.constantOutbreak1 = compute_information_criteria(FOIfit.constantOutbreak2)

seroprevalence.fit(FOIfit.constantOutbreak2, YLIM = 0.8) # plot model

FOIfit.constantOutbreak2$fit
pairs(FOIfit.constantOutbreak2$fit, pars = c("T[1]", "constant", "rho"))

stanObject_constantOutbreak2 <- FOIfit.constantOutbreak2$fit

## Check for convergence of chains, neff, Rhat and posteriors of all parameters
# Want ESS to be >200 
# Want Rhat <1.1
stan_diag(stanObject_constantOutbreak2)
stan_trace(stanObject_constantOutbreak2,pars = "T[1]")
stan_trace(stanObject_constantOutbreak2,pars = "rho")
stan_trace(stanObject_constantOutbreak2,pars = "constant")
stan_dens(stanObject_constantOutbreak2, pars = c("T[1]","rho","constant"))

launch_shinystan(stanObject_constantOutbreak2)

# save model
saveRDS(FOIfit.constantOutbreak2, "stan_constantOutbreak1_reverse_5_years")

FOIfit.constantOutbreak2 <- readRDS("savedModels/stan_constantOutbreak1_reverse_5_years")

## Extract log liklihood to calculate WAIC
logLikelihood_constantOutbreak2 <- extractLogLik(FOIfit.constantOutbreak2)

loo::waic(logLikelihood_constantOutbreak2)

#####################################################################
# 5. No constant FOI but with 1 outbreak in last 10 years
#####################################################################

outbreakModel = FOImodel( type='outbreak', ## specify type of model
                          K=1, ## number of outbreaks
                          seroreversion = 1, ## seroconversion is TRUE
                          priorT1 = 0, ## lower bound for outbreak time (time since sampling survey) (uniform distribution)
                          priorT2 = 5, ## upper bound for outbreak time (time since sampling survey) (uniform distribution)
                          priorRho1 = 0, ## lower bound for waning estimate (uniform distribution)
                          priorRho2 = 10) ## upper bound for waning estimate (uniform distribution)

FOIfit.outbreak = fit( data = seroModelData,  
                       model = outbreakModel, 
                       chains=4, 
                       iter = 20000)

parameters_credible_intervals(FOIfit.outbreak)
seroprevalence.fit(FOIfit.outbreak, YLIM = 0.8)

compute_information_criteria(FOIfit.outbreak)
FOIfit.outbreak$fit

pairs(FOIfit.outbreak$fit, pars = c("T[1]", "rho"))

stanObject_outbreak <- FOIfit.outbreak$fit

## Check for convergence of chains, neff, Rhat and posteriors of all parameters
# Want ESS to be >200 
# Want Rhat <1.1
stan_diag(stanObject_outbreak)
stan_trace(stanObject_outbreak,pars = "T[1]")
stan_trace(stanObject_outbreak,pars = "rho")
stan_dens(stanObject_outbreak, pars = c("T[1]","rho"))

launch_shinystan(stanObject_outbreak)

saveRDS(FOIfit.outbreak, "stan_outbreak1_reverse_10_years")

FOIfit.outbreak <- readRDS("savedModels/stan_outbreak1_reverse_10_years")

## Extract log liklihood to calculate WAIC
logLikelihood_outbreak <- extractLogLik(FOIfit.outbreak)

loo::waic(logLikelihood_outbreak)
