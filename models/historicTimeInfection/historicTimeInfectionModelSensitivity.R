#' Reconstructing time of historic infections
#' E Rees
#' Figure 4 (methods Figure 2)
#' Comparing population level - comparing seroprevalence from 2013 study with seroprevalence in 2012 outbreak
#' The 2012 distribution is the expected distribution of a recent outbreak
#' So for each individual in 2013 data, decide the probability with which they would be in each category
#' Sensitivity analysis - shifting the profiles from the 2012 dataset to match the geometric mean in Lupidi et al.
#' Due to data sharing constraints, the real data cannot be shared, instead 
#' the results have been created using DUMMY data

library(tidyverse)
library(gridExtra)
library(grid)
library(janitor)
library(viridis)
library(cowplot)
source("models/historicTimeInfection/historicTimeInfectionFunctions.R")

########################################################################
### Creating DUMMY data
########################################################################


## 1. Summary dataset of the 2013 seroprevalence MAT results (n=417)
## For individuals which are positive for more than one titre the highest titre is used

## Create dummy data 
## Frequency by titre 
titreLevels <- c(50,100,200,400,800,1600,3200) ## Different titre levels
nSero <- 417 ## Total number of positive individuals

# Randomly generate numbers of individuals positvie for each titre level
m <- matrix(runif(length(titreLevels),0,1), ncol=length(titreLevels))
m <- sweep(m, 1, rowSums(m), FUN="/")
m <- round(m *nSero)
m

# Create a dummy dataset of serovar titre, number positive, and frequency
resultsFreq2013 <- as_tibble(titreLevels) %>%
  mutate(n = m[1,]) %>%
  rename(seroMax = value) %>%
  mutate(freq = n/nSero) %>%
  mutate(seroMax = as.factor(seroMax))

## 2. 2012 Fiji suspected acute cases (n = 66)

nAcute = 66

# Randomly generate numbers of individuals positive for each titre level
m <- matrix(runif(length(titreLevels),0,1), ncol=length(titreLevels))
m <- sweep(m, 1, rowSums(m), FUN="/")
m <- round(m *nAcute)
m

# Create a dummy dataset of serovar titre, number positive, and frequency
resultsFreq2012 <- as_tibble(titreLevels) %>%
  mutate(n = m[1,]) %>%
  rename(seroMax = value) %>%
  mutate(freq = n/nAcute) %>%
  mutate(seroMax = as.factor(seroMax))

########################################################################

## Create shifted distributions, based on the geometric mean  
seroTitres <- c(50,100,200,400,800,1600,3200,6400)

## Create shifted distributions

seroN <- c(0,9,21,11,14,8,2,1)
resultsFreq2012_s1 <- data.frame(as.factor(seroTitres),seroN)
resultsFreq2012_s1 <- tbl_df(resultsFreq2012_s1) %>% 
  mutate(freq = (seroN / sum(resultsFreq2012_s1$seroN))) %>% 
  rename(seroMax =as.factor.seroTitres., n = seroN)

seroN <- c(0,0,9,21,11,14,8,3)
resultsFreq2012_s2 <- data.frame(as.factor(seroTitres),seroN)
resultsFreq2012_s2 <- tbl_df(resultsFreq2012_s2) %>% 
  mutate(freq = (seroN / sum(resultsFreq2012_s2$seroN))) %>% 
  rename(seroMax =as.factor.seroTitres., n = seroN)

seroN <- c(0,0,0,9,21,11,14,11)
resultsFreq2012_s3 <- data.frame(as.factor(seroTitres),seroN)
resultsFreq2012_s3 <- tbl_df(resultsFreq2012_s3) %>% 
  mutate(freq = (seroN / sum(resultsFreq2012_s3$seroN))) %>% 
  rename(seroMax =as.factor.seroTitres., n = seroN)

######################################################################## 
## Shifting distribution by 1 titre level
######################################################################## 

resultsFreq2012_s1 <- processAndShiftData(resultsFreq2012_s1)

newTallyDistribution <- process2013DatSA(resultsFreq2012_s1,resultsFreq2013)


newDistributionPlot <- ggplot(newTallyDistribution, aes(x=years,y=yearsSum)) +
  geom_bar(stat="identity", fill = "#31446BFF") +
  scale_x_discrete(breaks = c(0,0.66,1.32,1.98,2.64,3.3,3.96), labels = c("Q4 2013","Q1 2013","Q2 2012","Q4 2011","Q1 2011","Q2 2010", "Q4 2009"), name = "Likely time of infection") +
  scale_y_continuous(name = "Number of individuals") +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) 

######################################################################## 
## Shifting distribution by 2 titre levels
######################################################################## 

resultsFreq2012_s2 <- processAndShiftData(resultsFreq2012_s2)

newTallyDistribution <- process2013DatSA(resultsFreq2012_s2,resultsFreq2013)

newDistributionPlot2 <- ggplot(newTallyDistribution, aes(x=years,y=yearsSum)) +
  geom_bar(stat="identity", fill = "#31446BFF") +
  scale_x_discrete(breaks = c(0,0.66,1.32,1.98,2.64,3.3,3.96), labels = c("Q4 2013","Q1 2013","Q2 2012","Q4 2011","Q1 2011","Q2 2010", "Q4 2009"), name = "Likely time of infection") +
  scale_y_continuous(name = "Number of individuals") +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) 

######################################################################## 
## Shifting distribution by 3 titre levels
######################################################################## 

resultsFreq2012_s3 <- processAndShiftData(resultsFreq2012_s3)

newTallyDistribution <- process2013DatSA(resultsFreq2012_s3,resultsFreq2013)

newDistributionPlot3 <- ggplot(newTallyDistribution, aes(x=years,y=yearsSum)) +
  geom_bar(stat="identity", fill = "#31446BFF") +
  scale_x_discrete(breaks = c(0,0.66,1.32,1.98,2.64,3.3,3.96), labels = c("Q4 2013","Q1 2013","Q2 2012","Q4 2011","Q1 2011","Q2 2010", "Q4 2009"), name = "Likely time of infection") +
  scale_y_continuous(name = "Number of individuals") +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )

### Create a combined plot
## This will look different to the plots in the paper, since it was using dummy data

plot_grid(newDistributionPlot, newDistributionPlot2, newDistributionPlot3, labels = c('A', 'B',"C"),nrow = 1)
