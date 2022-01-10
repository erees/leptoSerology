#' Reconstructing time of historic infections
#' E Rees
#' Figure 4 (methods Figure 2)
#' Comparing population level - comparing seroprevalence from 2013 study with seroprevalence in 2012 outbreak
#' The 2012 distribution is the expected distribution of a recent outbreak
#' So for each individual in 2013 data, decide the probability with which they would be in each category
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

## 2. Individual level data assuming infections are independent events (n=520), by serovar
## individuals may appear more than once

## Create dummy data
nAll <- 520 # Sample size
serovars <- c("seroPohn","seroCopen","seroCanic","seroBallum","seroAustral","seroHad") ## Difference serovars

matTitres2013All <- as_tibble(seq(1:nAll)) %>%
  mutate(seroVal = sample(titreLevels,nAll,replace = TRUE)) %>% ## randomly sample from the titre levels
  mutate(serotype = sample(serovars,nAll,replace = TRUE)) %>% ## randomly sample from the serovats
  mutate(seroVal = as.factor(seroVal)) %>%
  rename(id = value) 


## 3. Individual level data using only the highest titre per individual (n=417), by serovar

## Create dummy data
matTitres2013 <- as_tibble(seq(1:nSero)) %>% ## id column
  mutate(seroMax = sample(titreLevels,nSero,replace = TRUE)) %>% ## randomly sample from the titre levels
  mutate(serotype = sample(serovars,nSero,replace = TRUE)) %>% ## randomly sample from the serovats
  mutate(seroMax = as.factor(seroMax)) %>%
  rename(id = value) 

## 4. 2012 Fiji suspected acute cases (n = 66)

nAcute = 66

# Randomly generate numbers of individuals positvie for each titre level
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
## Process 2012 data
######################################################################## 

resultsFreq2012 <- resultsFreq2012 %>%
  mutate(seroMax = as.character(seroMax)) %>%
  mutate(seroMax = as.numeric(seroMax)) %>%
  arrange(seroMax)

## Create a dataset where each new column sequentially shifts the column values down one row (increasing number of NA's at the beginning of each colum)
# A new proportion is then calculated for each row
## This provides an expected distribution of possible initial titre values, based on different MAT values from 2013
## i.e. If in 2013 an individual has a titre level of 1:400, we are interested in knowing the likely distribution of 1:400,1:800,1:1600 and 1:3200.

for(j in 1:6){
  resultsFreq2012[,ncol(resultsFreq2012) + 1] <- shiftDownAndFreq(resultsFreq2012$n, j)
  colnames(resultsFreq2012)[ncol(resultsFreq2012)] <- paste0("freq_", j)
}

resultsFreq2012

######################################################################## 
## Process 2013 data
######################################################################## 

titreDropYear <- 0.66 # (based on the estimated rate that individuals drop one antibody dilution level - 7.9 months)

## dataset with duplicates, assuming infections are independent events
mat2013BySerovarFull <- matTitres2013All %>%
  group_by(serotype) %>%
  dplyr::count(seroVal, sort = TRUE)

## dataset with no duplicates, assuming one one infecting titre per individual
mat2013BySerovarDeduped <- matTitres2013 %>%
  group_by(serotype) %>%
  dplyr::count(seroMax, sort = TRUE) %>%
  rename(seroVal = seroMax)

serovarNames <- c("seroPohn","seroCopen","seroAustral","seroCanic","seroBallum","seroHad")

## Function that wraps around the code so can repeat with full and de-duplicated data

produceSummaryDat <- function(mat2013BySerovar){
  
  ## Compare the MAT value in 2013 to the possible expected distribution of titres in 2012
  ## Repeat this for all serovars 
  
  seroPohnProcessed <- processSero2013("seroPohn",mat2013BySerovar)
  colnames(seroPohnProcessed) <- c("50","100","200","400","800","1600","3200","years")
  seroPohnProcessed <- tbl_df(seroPohnProcessed)
  newTallyDistributionPohn <- createTallyDat(seroPohnProcessed,"Pohnpei")
  
  seroCopenProcessed <- processSero2013("seroCopen",mat2013BySerovar)
  colnames(seroCopenProcessed) <- c("50","100","200","400","800","1600","3200","years")
  seroCopenProcessed <- tbl_df(seroCopenProcessed)
  newTallyDistributionCopen <- createTallyDat(seroCopenProcessed,"Copenhageni")
  
  seroAustralProcessed <- processSero2013("seroAustral",mat2013BySerovar)
  colnames(seroAustralProcessed) <- c("50","100","200","400","800","1600","3200","years")
  seroAustralProcessed <- tbl_df(seroAustralProcessed)
  newTallyDistributionAustral <- createTallyDat(seroAustralProcessed,"Australis")
  
  seroCanicProcessed <- processSero2013("seroCanic",mat2013BySerovar)
  colnames(seroCanicProcessed) <- c("50","100","200","400","800","1600","3200","years")
  seroCanicProcessed <- tbl_df(seroCanicProcessed)
  newTallyDistributionCanic <- createTallyDat(seroCanicProcessed,"Canicola")
  
  seroBallumProcessed <- processSero2013("seroBallum",mat2013BySerovar)
  colnames(seroBallumProcessed) <- c("50","100","200","400","800","1600","3200","years")
  seroBallumProcessed <- tbl_df(seroBallumProcessed)
  newTallyDistributionBallum <- createTallyDat(seroBallumProcessed,"Ballum")
  
  seroHadProcessed <- processSero2013("seroHad",mat2013BySerovar)
  colnames(seroHadProcessed) <- c("50","100","200","400","800","1600","3200","years")
  seroHadProcessed <- tbl_df(seroHadProcessed)
  newTallyDistributionHad <- createTallyDat(seroHadProcessed,"Hardjo")
  
  serovarResults <- newTallyDistributionPohn %>%
    rbind(newTallyDistributionAustral) %>%
    rbind(newTallyDistributionCanic) %>%
    rbind(newTallyDistributionCopen) %>%
    rbind(newTallyDistributionBallum) %>%
    rbind(newTallyDistributionHad)
  
  ## Can either output frequencies or absolute numbers
  serovarResultsTest <- serovarResults %>%
    spread(serovar, yearsSum)
  
  serovarResultsFreq <- serovarResultsTest %>%
    janitor::adorn_percentages("col") 
  
  serovarResultsFreq <- serovarResultsFreq %>%
    gather("serovar","n",2:7)
  
  serovarResultsFreq <- serovarResultsFreq %>%
    mutate(years = as.factor(years))
  
  serovarResultsTest <- serovarResultsTest %>%
    gather("serovar","n",2:7)
  
  serovarResultsTest <- serovarResultsTest %>%
    mutate(years = as.factor(years))
  
  serovarResultsTest
}

serovarResultsTestFull <- produceSummaryDat(mat2013BySerovarFull)

serovarResultsTestDeduped <- produceSummaryDat(mat2013BySerovarDeduped)

createPlot <- function(serovarResultsTest){
  ggplot(serovarResultsTest, aes(x=years,y=n, fill = serovar)) +
    geom_bar(stat="identity") +
    scale_x_discrete(breaks = c(0,0.66,1.32,1.98,2.64,3.3,3.96), labels = c("03-2013 -\n11-2013\n","07-2012 -\n02-2013\n","11-2011 -\n06-2012\n","03-2011 -\n10-2011\n","07-2010 -\n02-2011\n", "11-2009 -\n07-2010\n", "03-2009 -\n10-2009\n"), name = "Likely time of infection") +
    scale_fill_viridis(discrete = TRUE, option="cividis",direction = -1, "Serovars") +
    scale_y_continuous(name = "Number of individuals",limits = c(0,155)) +
    theme(legend.title = element_blank(), 
          # legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          axis.text.x=element_text(size=7)
    ) 
}

plotA <- createPlot(serovarResultsTestFull)
plotB <- createPlot(serovarResultsTestDeduped)

## This will look different to the plots in the paper, since it was using dummy data
plot_grid(plotA, plotB,labels = "AUTO")

