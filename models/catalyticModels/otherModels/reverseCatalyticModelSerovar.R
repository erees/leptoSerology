############################################################################
## Reverse catalytic model by division
## Analyse data by serovar - allowing FOI to vary by serovar, waning held constant
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
require(MCMCvis)
require(loo)

############################################################################
# Read in and process full data
# Consists of line level data with serostatus, age, and serovar
# Data not available for sharing
#
# Split the data by serovar (i.e. end up with 5 different datasets, 
# with the serostatus by serovar)
# Put them into a long dataset, since estimating waning across all serovars
############################################################################



############################################################################
# Model reverse catalytic model by serovar
# allows foi to vary by serovar, but keeps waning the same accross all. 
############################################################################

jcode <- "model{
for (i in 1:2150){ 
n.pos[i] ~ dbern(seropos_est[i]) 
seropos_est[i] = (lambdaPohnpei / (lambdaPohnpei + delta)) * (1-exp(-(lambdaPohnpei+delta)*age[i])) #reverse catalytic model
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)
}
for (i in 2151:4300){ 
n.pos[i] ~ dbern(seropos_est[i]) 
seropos_est[i] = (lambdaCopen / (lambdaCopen + delta)) * (1-exp(-(lambdaCopen+delta)*age[i])) 
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)
}
for (i in 4301:6450){ 
n.pos[i] ~ dbern(seropos_est[i]) 
seropos_est[i] = (lambdaAustral / (lambdaAustral + delta)) * (1-exp(-(lambdaAustral+delta)*age[i])) 
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)
}
for (i in 6451:8600){ 
n.pos[i] ~ dbern(seropos_est[i]) 
seropos_est[i] = (lambdaBallum / (lambdaBallum + delta)) * (1-exp(-(lambdaBallum+delta)*age[i])) 
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)
}
for (i in 8601:10750){ 
n.pos[i] ~ dbern(seropos_est[i]) 
seropos_est[i] = (lambdaCanic / (lambdaCanic + delta)) * (1-exp(-(lambdaCanic+delta)*age[i])) 
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)
}

# priors
lambdaPohnpei ~ dunif(0,0.1) 
lambdaCopen ~ dunif(0,0.1) 
lambdaAustral ~ dunif(0,0.1) 
lambdaBallum ~ dunif(0,0.1) 
lambdaCanic ~ dunif(0,0.1) 
delta ~ dunif(0,10) 
}"

############################################################################
# Run model
############################################################################

mcmc.length=10000

jdat <- list(n.pos= serovarDf$value,
             age=serovarDf$age)
jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=5000)
update(jmod)
jpos = coda.samples(jmod, c("lambdaPohnpei","lambdaCopen","lambdaAustral", "lambdaBallum", "lambdaCanic","delta","loglik"), n.iter=mcmc.length)

jpos1 <- jpos[,1:6,drop=TRUE]
plot(jpos1) # check convergence
summary(jpos1)
MCMCsummary(jpos1, round = 3)


## check convergence
plot(jpos)
summary(jpos)

autocorr.plot(jpos)

jpos <- window(jpos, start=12001)
summary(jpos)

dic.samples(jmod, n.iter = mcmc.length) ##calculate DIC

mcmcMatrix <- as.matrix(jpos)

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

## extract log liklihood
logLik <- mcmcMatrix[,7:2156]

## Caclulate WAIC and LOO from loo package
waic <- waic(logLik)
waic
loo <- loo(logLik)
loo
plot(loo,label_points = TRUE)