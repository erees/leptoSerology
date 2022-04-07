# leptoSerology

This repository contains the data and code used to produce the results presented in "Estimating the duration of antibody positivity and likely time of Leptospira infection using data from a cross-sectional serological study in Fiji".

Due to data sharing constraints, the full dataset could not be shared. Instead, aggregated data by 5-year age groups is provided in the data folder (`seroDat_grouped5.csv`). 

The main models from the paper (catalytic and reverse catalytic model) have been modified to use this data, and can be found here:
* `models\catalyticModels\catalyticModel.R` 
* `models\catalyticModels\reverseCatalyticModel.R` 

The results from these models are very similar to those presented in the paper, although there will be slight differences in the estimates. The code to recreate figure 3 can be found here:
* `models\catalyticModels\figure3.R`

The catalytic models by sex, division, serovar, and time varying FOI use the full dataset, and therefor cannot be re-created. For completeness, the code used to produce these resuls (without the data) are included within `models\catalyticModels\otherModels`.

The simulation recovery study (used to create supplementary figures 4 and 5 and supplementary table 3) can be found here:
`models\catalyticModels\otherModels\simulationRecoveryStudy.R`

The code to estimate the historic time of infection is within `models\historicTimeInfection` folder. 

Firstly, the antibody titre drop was estimated from Lupidi _et al_. The data extracted from the paper is in `data\Lupidi_full_data.csv`. The code to estimate the antibody titre drop can be found here:
* `models\historicTimeInfection\estimatingAntibodyDrop.R`.

To estimate the histric time of infection from the seroprevalence survey, individual MAT titres were used. These were not available to share, so instead, dummy data has been created to demonstrate how this was estimated. The code for this can befound here:
* `models\historicTimeInfection\historicTimeInfectionModel.R`
* `models\historicTimeInfection\historicTimeInfectionSensitivity.R` (Sensitivity analysis)
