############################################################################
## Reverse catalytic model by sex
## Analyse data by sex - allowing FOI to vary by sex, waning held constant
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
# Consists of line level data with serostatus, age, and sex
## Data not available for sharing
############################################################################


############################################################################
### define model
### Reverse catalytic model by sex
### FOI estimated by division, waning held across sexes
############################################################################

jcode <- "model{ 
#Females
for (i in 1:1160){ 
n.pos[i] ~ dbern(seropos_est[i])
seropos_est[i] = (lambdaFem / (lambdaFem + delta)) * (1-exp(-(lambdaFem+delta)*age[i])) #reverse catalytic model
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)
}
# Males
for (i in 1161:length(n.pos)){ 
n.pos[i] ~ dbern(seropos_est[i])
seropos_est[i] = (lambdaMal / (lambdaMal + delta)) * (1-exp(-(lambdaMal+delta)*age[i])) #reverse catalytic model
loglik[i] <- logdensity.bin(n.pos[i],seropos_est[i],1)
}

# Was getting error so had to make prior smaller (dunif(0,1) was too wide)
lambdaFem ~ dunif(0,0.5) #uninformative prior
lambdaMal ~ dunif(0,0.5) #uninformative prior
delta ~ dunif(0,10) #uniformative prior 
}"


################################################################################
# Run the model
################################################################################

mcmc.length=60000

jdat <- list(n.pos= seroDataBin$seropositive,
             age=seroDataBin$age)

jmod = jags.model(textConnection(jcode), 
                  data=jdat, 
                  n.chains=4, 
                  n.adapt=10000)
update(jmod)

## Includes likelihood for each data point 
jpos = coda.samples(jmod, c("lambdaFem", "lambdaMal","delta", "loglik"), n.iter=mcmc.length)

## Extract only the first three columns - includes only the parameters
jpos1 <- jpos[,1:3,drop=TRUE]
plot(jpos1) # check convergence
summary(jpos1)
jpos1 <- window(jpos1, start=50001)

MCMCsummary(jpos1, round = 2) ## Check ESS and Rhat

#convert mcmc.list to a matrix
mcmcMatrix <- as.matrix(jpos1)


# Plotting posterior distributions of all parameters
mcmcDF <- as_tibble(mcmcMatrix)
mcmcDF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

# calculate the DIC 
dic.samples(jmod, n.iter = mcmc.length)

#convert mcmc.list to a matrix
mcmcMatrixlh <- as.matrix(jpos)

## extract log liklihood
logLik <- mcmcMatrixlh[,4:2153]

## Caclulate WAIC and LOO from loo package
waic <- waic(logLik)
waic
loo <- loo(logLik)
loo
plot(loo,label_points = TRUE)