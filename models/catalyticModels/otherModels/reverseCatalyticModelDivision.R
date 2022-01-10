############################################################################
## Reverse catalytic model by division
## Analyse data by sex - allowing FOI to vary by division, waning held constant
## E Rees
##
## FOR REFERENCE
## Due to data sharing constraints it is not possible to share the underlying data for this code
## This model is provided for completeness
############################################################################

require(tidyverse)
require(rjags)
require(binom)
require(varhandle)
require(loo)
require(MCMCvis)

############################################################################
# Read in and process full data
# Consists of line level data with serostatus, age, and division
## Data not available for sharing
############################################################################

############################################################################
### define model
### Reverse catalytic model by division
### FOI estimated by division, waning held across all divisions
############################################################################

jcode <- "model{
#central
for (i in 1:662){ 
n.pos[i] ~ dbern(seropos_est[i])
seropos_est[i] = (lambdaCentral / (lambdaCentral + delta)) * (1-exp(-(lambdaCentral+delta)*age[i])) #reverse catalytic model
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)

}
# north
for (i in 663:1177){ 
n.pos[i] ~ dbern(seropos_est[i])
seropos_est[i] = (lambdaNorth / (lambdaNorth + delta)) * (1-exp(-(lambdaNorth+delta)*age[i])) #reverse catalytic model
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)

}
# west
for (i in 1178:length(n.pos)){ 
n.pos[i] ~ dbern(seropos_est[i])
seropos_est[i] = (lambdaWest / (lambdaWest + delta)) * (1-exp(-(lambdaWest+delta)*age[i])) #reverse catalytic model
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)

}

# priors
lambdaCentral ~ dunif(0,0.5) #uninformative prior
lambdaNorth ~ dunif(0,0.5) #uninformative prior
lambdaWest ~ dunif(0,0.5) #uninformative prior

delta ~ dunif(0,10) #uniformative prior 

}"

################################################################################
# Running the model
################################################################################




mcmc.length=60000

jdat <- list(n.pos= seroDataBin$seropositive,
             age=seroDataBin$age)
jmod = jags.model(textConnection(jcode), 
                  data=jdat, 
                  n.chains=4, 
                  n.adapt=20000)
update(jmod)
jpos = coda.samples(jmod, c("lambdaCentral", 
                            "lambdaNorth", 
                            "lambdaWest",
                            "delta",
                            "loglik"), 
                    n.iter=mcmc.length)

jpos1 <- jpos[,1:4,drop=TRUE]
plot(jpos1) # check convergence
summary(jpos1)
MCMCsummary(jpos1, round = 3)

plot(jpos) # check convergence
summary(jpos)
jpos1 <- window(jpos1, start=60001)

MCMCsummary(jpos, round = 2) ## Check ESS and Rhat

jpos <- window(jpos, start=60001)

#convert mcmc.list to a matrix
mcmcMatrix <- as.matrix(jpos)

# Plotting posterior distributions of all parameters
mcmcDF <- as_tibble(mcmcMatrix)
mcmcDF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

# calculate the DIC 
dic.samples(jmod, n.iter = mcmc.length)

## extract log liklihood
logLik <- mcmcMatrix[,5:2154]

## Caclulate WAIC and LOO from loo package
waic <- waic(logLik)
waic
loo <- loo(logLik)
loo
plot(loo,label_points = TRUE)