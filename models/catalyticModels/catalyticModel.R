############################################################################
##  Catalytic model 
## E Rees
## Due to data sharing constrains, this code used data binned into 5 year age 
## categories rather than the individual data
## Therefore, the results are slightly different to those that appear in the paper
############################################################################

## Read in packages
require(tidyverse)
require(rjags)
require(binom)
require(varhandle)
require(MCMCvis)
require(loo)

################################################################################
### Read in data and source functions
################################################################################

# read in model functions
source("models/catalyticModels/catalyticModelfunctions.R")

# Reading in 5-year age grouped data
seroDatGrouped <- read_csv("data/seroDat_grouped5.csv")

fiveYearAgeTotals <- seroDatGrouped$total

######################################
### define model catalytic model
######################################

#define model
jcode <- "model{ 
for (i in 1:length(n)){

n.pos[i] ~ dbin(seropos_est[i],n[i]) 

#catalytic model
seropos_est[i] = 1-exp(-lambda*age[i]) 

#calculate likelihood (used for WAIC)
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],n[i]) 

}
## Priors
lambda ~ dunif(0,0.5) 
}"


################################################################################
# Running the model
################################################################################

## Number of model itterations
mcmc.length=10000

## Specify my data
jdat <- list(n.pos = seroDatGrouped$seropositive,
             age = seroDatGrouped$midpoint,
             n=seroDatGrouped$total)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=2000)
update(jmod)
jpos = coda.samples(jmod, c("lambda","loglik"), n.iter=mcmc.length)

plot(jpos) # check convergence

summary(jpos)
MCMCsummary(jpos, round = 2) ## Check ESS and Rhat

#convert mcmc.list to a matrix
mcmcMatrix <- as.matrix(jpos)

# Plotting posterior distributions of all parameters
mcmcDF <- as_tibble(mcmcMatrix)
mcmcDF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

# Get point estimates
lambdaPointEst <- mcmcMatrix[,"lambda"] %>% quantile(probs=c(.5,.025,.975))

# calculate DIC 
dic.samples(jmod, n.iter = mcmc.length)

## extract just the log liklihood
logLik <- mcmcMatrix[,2:16]

## Caclulate WAIC and LOO from loo package
waic <- waic(logLik)
waic
loo <- loo(logLik)
loo
plot(loo,label_points = TRUE)


################################################################################
## Get data in correct format for plotting
## Sample from mcmc chains for credible intervals
## Add binomial sampling uncertainty
################################################################################

midpoints <- seroDatGrouped$midpoint
ager=0:80

## 1. Sample from mcmc chain to get 95% credible intervals (model uncertainty)
numSamples = 1000
outDf <- matrix(NA, nrow=numSamples, ncol=1)

for (i in 1:numSamples ) {
  randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
  outDf[i,] <- mcmcMatrix[randomNumber,][[1]]
}

lambdaEst <- outDf %>% quantile(probs=c(.5,.025,.975))
lambdaEst %>% round(3) %>% print

## Create a df with model uncertainty
df_mod=data.frame(midpoint=ager, 
                  mean=1-exp(-lambdaEst[1]*ager),
                  lower=1-exp(-lambdaEst[2]*ager),
                  upper=1-exp(-lambdaEst[3]*ager))


## 2. Binomial sample uncertainty - accounts for the sample size of the underlying data
randomlySampleMcmcChain <- mcmcRandomSamplerCat(1000,mcmcMatrix,midpoints,fiveYearAgeTotals)
ageQuantilesSamplingUncertainty <- ageQuantiles(randomlySampleMcmcChain)

## Create a df with sample uncertainty
df_sampling = data.frame(
  midpoint = seroDatGrouped$midpoint,
  mean = 1 - exp(-seroDatGrouped$midpoint*(lambdaEst[1])),
  upper = ageQuantilesSamplingUncertainty[,3],
  lower = ageQuantilesSamplingUncertainty[,2]
)

## Save the dataframes for plotting (along with reverse catalytic model)
saveRDS(df_sampling,"samplingUncertaintyCat_5.rds")
saveRDS(df_mod, "modelUncertaintyCat_5.rds" )

