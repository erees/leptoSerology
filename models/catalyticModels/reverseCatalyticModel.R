############################################################################
## Reverse catalytic model 
## E Rees
## Due to data sharing constraints, this code used data binned into 5 year age 
## categories rather than the individual data
## Therefore, the results are slightly different to those that appear in the paper
############################################################################

## Load in packages
require(tidyverse)
require(rjags)
require(binom)
require(varhandle)
require(loo)
require(MCMCvis)

################################################################################
### Read in data and source functions
################################################################################

# read in model functions
source("models/catalyticModels/catalyticModelfunctions.R")

# Reading in 5-year age grouped data
seroDatGrouped <- read_csv("data/seroDat_grouped5.csv")

fiveYearAgeTotals <- seroDatGrouped$total


###################
### define model
###################

jcode <- "model{ 
for (i in 1:length(n)){

n.pos[i] ~ dbin(seropos_est[i],n[i]) 

#reverse catalytic model
seropos_est[i] = (lambda / (lambda + delta)) * (1-exp(-(lambda+delta)*age[i])) 

## #calculate likelihood (used for WAIC)
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],n[i]) 
}
# Define priors
lambda ~ dunif(0,0.5)
delta ~ dunif(0,10)
}"

################################################################################
# Running the model
################################################################################

## Number of model itterations
mcmc.length=20000 

## Specify my data
jdat <- list(n.pos= seroDatGrouped$seropositive,
             age=seroDatGrouped$midpoint,
             n=seroDatGrouped$total)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=5000)
update(jmod)
jpos = coda.samples(jmod, c("lambda","delta","loglik"), n.iter=mcmc.length)

plot(jpos) ## Check convergence of all chains

MCMCsummary(jpos, round = 2) ## Check ESS and Rhat
summary(jpos)

#convert mcmc.list to a matrix
mcmcMatrix <- as.matrix(jpos)

# Plotting posterior distributions of all parameters
mcmcDF <- as_tibble(mcmcMatrix)
mcmcDF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

#remove burn in
jpos <- window(jpos, start=10000)
plot(jpos)
MCMCsummary(jpos, round = 2)
summary(jpos)

# calculate the DIC 
dic.samples(jmod, n.iter = mcmc.length)

## extract log liklihood
logLik <- mcmcMatrix[,3:2152]

## Caclulate WAIC and LOO from loo package
waic <- waic(logLik)
waic
loo <- loo(logLik)
loo
plot(loo,label_points = TRUE)

# Get point estimates
deltaPointEst <- mcmcMatrix[,"delta"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEst <- mcmcMatrix[,"lambda"] %>% quantile(probs=c(.5,.025,.975))

################################################################################
## Get data in correct format for plotting
## Sample from mcmc chains for credible intervals
## Add sampling uncertainty
################################################################################

ageVector = 0:80
numberOfSamples = 1000
outDf <- matrix(NA,nrow=numberOfSamples, ncol = length(ageVector))
midpoints <- seroDatGrouped$midpoint


## 1. Sample from mcmc chain to get 95% credible intervals (model uncertainty)

randomlySampleMcmcChainModelUncertainty <- mcmcModelRandomSampler(1000,mcmcMatrix,ageVector)
# for each column in the matrix get quantiles by age
ageQuantilesModelUncertainty <- ageQuantiles(randomlySampleMcmcChainModelUncertainty)

## Create a df with model uncertainty
df_upperLower = data.frame(
  midpoint = ageVector,
  mean = (lambdaPointEst[1] / (lambdaPointEst[1]+deltaPointEst[1])) *(1 - exp(-ageVector*(lambdaPointEst[1]+deltaPointEst[1]))),
  upper = ageQuantilesModelUncertainty[,3],
  lower = ageQuantilesModelUncertainty[,2]
)

## 2. Sample uncertainty - accounts for the sample size of the underlying data
randomlySampleMcmcChain <- mcmcRandomSampler(1000,mcmcMatrix,midpoints,fiveYearAgeTotals)
ageQuantilesSamplingUncertainty <- ageQuantiles(randomlySampleMcmcChain)

## Create a df with sample uncertainty
df_sampling = data.frame(
  midpoint = seroDatGrouped$midpoint,
  mean = (lambdaPointEst[1] / (lambdaPointEst[1]+deltaPointEst[1])) *(1 - exp(-seroDatGrouped$midpoint*(lambdaPointEst[1]+deltaPointEst[1]))),
  upper = ageQuantilesSamplingUncertainty[,3],
  lower = ageQuantilesSamplingUncertainty[,2]
)

## Save the dataframes for plotting
saveRDS(df_sampling,"samplingUncertaintyRevCat_5.rds")
saveRDS(df_upperLower, "modelUncertaintyRevCat_5.rds")



