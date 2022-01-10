############################################################################
##  Code used to create Fig 3 
## Catalytic and reverse catalytic model outputs
## E Rees
############################################################################

require(tidyverse)
require(cowplot)

## Read in dataframes for catalytic and reverse catalytic models (sampling and model uncertainty)
samplingUncertaintyCat <- readRDS("savedModelOutputs/samplingUncertaintyCat_5.rds")
modelUncertaintyCat <- readRDS("savedModelOutputs/modelUncertaintyCat_5.rds" )

samplingUncertaintyRevCat <- readRDS("savedModelOutputs/samplingUncertaintyRevCat_5.rds")
modelUncertaintyRevCat <- readRDS("savedModelOutputs/modelUncertaintyRevCat_5.rds" )

## Read in 5 year age grouped data
seroDat <- readRDS("data/seroDatGrouped.rds") ## seroprevalence binned into 5 year age groups (used for plotting)

## Reverse catalytic plot
revCatPlot <- ggplot(modelUncertaintyRevCat, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.3, fill = "#457b9d")+
  geom_line()+
  geom_point(data=seroDat)+
  geom_linerange(data=seroDat) +
  geom_ribbon(data=samplingUncertaintyRevCat, alpha=0.3, fill = "#457b9d")+
  scale_y_continuous(breaks=seq(0,0.6,by=0.1), limits = c(0,0.55))+
  xlab("Age (years)") + ylab("Proportion seropositive") +
  scale_x_continuous(breaks=seq(0,100,by=10)) +
  theme_minimal()

## Catalytic plot
catPlot <- ggplot(modelUncertaintyCat, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.3, fill = "#9e2a2b")+
  geom_line()+
  geom_point(data=seroDat)+
  geom_linerange(data=seroDat) +
  geom_ribbon(data=samplingUncertaintyCat, alpha=0.3,fill = "#9e2a2b")+
  scale_y_continuous(breaks=seq(0,0.6,by=0.1), limits = c(0,0.55))+
  xlab("Age (years)") + ylab("Proportion seropositive") +
  scale_x_continuous(breaks=seq(0,80,by=10)) +
  theme_minimal()

## Create panel plot
plot_grid(catPlot, revCatPlot, labels = c('A', 'B'))
