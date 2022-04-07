############################################################################
## Two scenarios, one with high FOI and one with low FOI
## Waning the same in both settings
## Simulation study to see whether FOI and waning can be recovered

## Possible to recover FOI and waning more accurately in the high FOI setting 
## but not the low FOI setting

## E Rees
## 06/04/2022
############################################################################


## Read in packages
library(tidyverse)
library(binom)
library(varhandle)
library(rjags)
library(MCMCvis)
library(cowplot)


## set ages
ager=0:80

## set parameter values
deltaSample <- 0.1
lambdaHigh <- 0.05
lambdaLow <- 0.005

## reverse catalytic model for high and low FOI setting
highFOI <- (lambdaHigh / (lambdaHigh+deltaSample)) * (1 - exp(-ager*(lambdaHigh+deltaSample)))
lowFOI <- (lambdaLow / (lambdaLow+deltaSample)) * (1 - exp(-ager*(lambdaLow+deltaSample)))

## Add these to a df for plotting
outDf <- data.frame(age = ager,seroHigh = highFOI, seroLow = lowFOI)
outDf <- as_tibble(outDf)


### From each model, sample 50 times for every age (0-80 years)
total = 50 

sampleBinom <- function(prob){
  sum(rbinom(n = total,size = 1,prob = prob))
}

outDfHigh <- outDf %>%
  mutate(seroSample = sapply(X = outDf$seroHigh, sampleBinom)) %>%
  mutate(total = total)

outDfLow <- outDf %>%
  mutate(seroSample = sapply(X = outDf$seroLow, sampleBinom)) %>%
  mutate(total = total)

## categorise sampled data into 5-year age groups
ageCuts <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,81)
ageLabels <- c(2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5)

seroDatHigh <- outDfHigh %>%
  mutate(ageCat = cut(age,ageCuts , labels=ageLabels)) %>%
  mutate(midpoint = unfactor(ageCat)) %>%
  group_by(midpoint) %>%
  summarise(total =sum(total), seropositive = sum(seroSample)) %>%
  slice(-17) %>%
  mutate(setting = "High")

seroDatLow <- outDfLow %>%
  mutate(ageCat = cut(age, ageCuts, labels=ageLabels)) %>%
  mutate(midpoint = unfactor(ageCat)) %>%
  group_by(midpoint) %>%
  summarise(total =sum(total), seropositive = sum(seroSample)) %>%
  slice(-17) %>%
  mutate(setting = "Low")

## Calculate binomial CI for each 5-year age group
seroDatHigh[,c("mean","lower","upper")] <- binom.confint(seroDatHigh$seropositive, seroDatHigh$total,method="exact")[,c("mean","lower","upper")]
seroDatLow[,c("mean","lower","upper")] <- binom.confint(seroDatLow$seropositive, seroDatLow$total,method="exact")[,c("mean","lower","upper")]

## Combine into one dataset for plotting
seroDat <- seroDatHigh %>%
  rbind(seroDatLow)

## Create plot (supplementary figure 4) showing the reverse catalytic model (line)
## along with the samples mean and CI for each 5-year age group
## For each setting (high and low FOI)
ggplot(outDf, aes(x=age, y=seroHigh)) +
  geom_line(color = "darkgoldenrod")+
  geom_line(aes(y=seroLow),color = "darkcyan")+
  geom_point(data=seroDat,aes(x=midpoint,y=mean, color = setting), position = position_dodge(width = 1))+
  geom_linerange(data=seroDat, aes(x = midpoint,y=mean, ymin = lower, ymax = upper, color=setting),position = position_dodge(width = 1)) +
  scale_colour_manual(name = "Setting", labels = c("High","Low"), values = c("darkgoldenrod","darkcyan")) +
  xlab("Age (years)") + ylab("Proportion seropositive") +
  scale_x_continuous(breaks=seq(0,80,by=10)) +
  theme_minimal()

## Now estimate the FOI and waning from the sampled seroprevalence estimates
## To see if we can recover the known parameters

#define model to do this
jcode <- "model{ 
for (i in 1:length(n.pos)){
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] = (lambda / (lambda + delta)) * (1-exp(-(lambda+delta)*age[i])) #reverse catalytic model

}

lambda ~ dgamma(1.1,3) #gamma prior (mean and rate)
delta ~ dnorm(1.1,0.3) #gamma prior
}"



############################################################################
##### High FOI setting
############################################################################

# Run model

mcmc.length=120000 # number of mcmc samples
jdat <- list(n.pos= seroDatHigh$seropositive,
             N=seroDatHigh$total,
             age=seroDatHigh$midpoint)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=20000)
update(jmod)
jpos = coda.samples(jmod, c("lambda","delta"), n.iter=mcmc.length,thin = 5)

# Gives you ESS and rhat - and parameter estimates
## Parameter estimates for supplemetary table 3
MCMCsummary(jpos, round = 3)

mcmcMatrix <- as.matrix(jpos)
mcmcDf <- as_tibble(mcmcMatrix)

p1 <- ggplot(mcmcDf, aes(x = delta)) +
  geom_density(alpha=0.50,size=0.3,fill="darkgoldenrod") + 
  theme_bw() +
  xlab("Waning Estimate") +
  ylab("Density") 


p2 <- ggplot(mcmcDf, aes(x = lambda)) +
  geom_density(alpha=0.50,size=0.3,fill="darkgoldenrod") + 
  theme_bw() +
  xlab("FOI Estimate") +
  ylab("Density") 


############################################################################
##### Low FOI setting
############################################################################

# Run model

mcmc.length=120000
jdat <- list(n.pos= seroDatLow$seropositive,
             N=seroDatLow$total,
             age=seroDatLow$midpoint)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=20000)
update(jmod)
jpos = coda.samples(jmod, c("lambda","delta"), n.iter=mcmc.length,thin=5)

# Gives you ESS and rhat
## Parameter estimates for supplemetary table 3
MCMCsummary(jpos, round = 3)

mcmcMatrix <- as.matrix(jpos)
mcmcDf <- as_tibble(mcmcMatrix)

p3 <- ggplot(mcmcDf, aes(x = delta)) +
  geom_density(alpha=0.50,size=0.3,fill="darkcyan") + 
  theme_bw() +
  xlab("Waning Estimate") +
  ylab("Density") 


p4 <- ggplot(mcmcDf, aes(x = lambda)) +
  geom_density(alpha=0.50,size=0.3,fill="darkcyan") + 
  theme_bw() +
  xlab("FOI Estimate") +
  ylab("Density") 


## Create supplementary Figure 5
plot_grid(p1,p2,p3,p4,ncol=2,labels="AUTO")

