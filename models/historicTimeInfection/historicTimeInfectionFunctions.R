
## Function which shifts down the column N times adding NA's, and then calculates a new proportion for the vector
shiftDownAndFreq <- function(x, n){
  shiftedVec <- c(rep(NA, n), tail(x, -n))
  newVec <- c()
  for(i in 1:length(shiftedVec)) {
    prop <- shiftedVec[i] / sum(shiftedVec,na.rm = TRUE)
    newVec[i] <- prop
  }
  newVec
}

## create a function to process the serovar data
processSero2013 <- function(serovarName,mat2013BySerovar){
  
  matSerovarSpec <- mat2013BySerovar %>% filter(serotype == serovarName)
  serovarMatrix <- matrix(,nrow = 7, ncol = 0) ## create an empty matrix with 7 rows
  
  if(50 %in% matSerovarSpec$seroVal) {
    matSerovarSpec50 <- matSerovarSpec[matSerovarSpec$seroVal == 50,]$n * resultsFreq2012$freq
    serovarMatrix <- cbind(serovarMatrix,matSerovarSpec50)
  }
  
  if(100 %in% matSerovarSpec$seroVal) {
    matSerovarSpec100 <- matSerovarSpec[matSerovarSpec$seroVal == 100,]$n * resultsFreq2012$freq_1
    matSerovarSpec100 <- c(matSerovarSpec100[2:7],NA)
    serovarMatrix <- cbind(serovarMatrix,matSerovarSpec100)
  }
  
  if(200 %in% matSerovarSpec$seroVal) {
    matSerovarSpec200 <- matSerovarSpec[matSerovarSpec$seroVal == 200,]$n * resultsFreq2012$freq_2
    matSerovarSpec200 <- c(matSerovarSpec200[3:7],NA,NA)
    serovarMatrix <- cbind(serovarMatrix,matSerovarSpec200)
  }
  
  if(400 %in% matSerovarSpec$seroVal) {
    matSerovarSpec400 <- matSerovarSpec[matSerovarSpec$seroVal == 400,]$n * resultsFreq2012$freq_3
    matSerovarSpec400 <- c(matSerovarSpec400[4:7],NA,NA,NA)
    serovarMatrix <- cbind(serovarMatrix,matSerovarSpec400)
  }
  
  if(800 %in% matSerovarSpec$seroVal) {
    matSerovarSpec800 <- matSerovarSpec[matSerovarSpec$seroVal == 800,]$n * resultsFreq2012$freq_4
    matSerovarSpec800 <- c(matSerovarSpec800[5:7],NA,NA,NA,NA)
    serovarMatrix <- cbind(serovarMatrix,matSerovarSpec800)
  }
  
  if(1600 %in% matSerovarSpec$seroVal) {
    matSerovarSpec1600 <- matSerovarSpec[matSerovarSpec$seroVal == 1600,]$n * resultsFreq2012$freq_5
    matSerovarSpec1600 <- c(matSerovarSpec1600[6:7],NA,NA,NA,NA,NA)
    serovarMatrix <- cbind(serovarMatrix,matSerovarSpec1600)
  }
  
  if(3200 %in% matSerovarSpec$seroVal) {
    matSerovarSpec3200 <- matSerovarSpec[matSerovarSpec$seroVal == 3200,]$n * resultsFreq2012$freq_6
    matSerovarSpec3200 <- c(matSerovarSpec3200[7],NA,NA,NA,NA,NA,NA)
    serovarMatrix <- cbind(serovarMatrix,matSerovarSpec3200)
  }
  
  serovarMatrix <- cbind(serovarMatrix,c(0,1*titreDropYear,2*titreDropYear,3*titreDropYear,4*titreDropYear,5*titreDropYear,6*titreDropYear))
  serovarMatrix
  
}

## Function which creates tally distribution
createTallyDat <- function(seroDatProcessed, seroName){
  
  newTallyDistribution <- seroDatProcessed %>%
    replace(., is.na(.), 0) %>%
    mutate(yearsSum = rowSums(.[1:(ncol(seroDatProcessed)-1)])) %>%
    mutate(years = as.factor(years))
  
  totalYearSum <- sum(newTallyDistribution$yearsSum)
  
  newTallyDistribution <- newTallyDistribution %>%
    mutate(freq = yearsSum/totalYearSum) %>%
    mutate(serovar = seroName) %>%
    select(years, yearsSum, serovar)
  
  newTallyDistribution
}

processAndShiftData <- function(resultsFreq2012){
  resultsFreq2012 <- resultsFreq2012 %>%
    mutate(seroMax = as.character(seroMax)) %>%
    mutate(seroMax = as.numeric(seroMax)) %>%
    arrange(seroMax)
  
  
  resultsFreq2012_1 <- tail(resultsFreq2012,-1) %>%
    mutate(freq_1 = n / sum(n)) 
  resultsFreq2012_1 <- c(NA,resultsFreq2012_1$freq_1)
  
  resultsFreq2012_2 <- tail(resultsFreq2012,-2) %>%
    mutate(freq_1 = n / sum(n)) 
  resultsFreq2012_2 <- c(NA,NA,resultsFreq2012_2$freq_1)
  
  resultsFreq2012_3 <- tail(resultsFreq2012,-3) %>%
    mutate(freq_1 = n / sum(n)) 
  resultsFreq2012_3 <- c(NA,NA,NA,resultsFreq2012_3$freq_1)
  
  resultsFreq2012_4 <- tail(resultsFreq2012,-4) %>%
    mutate(freq_1 = n / sum(n)) 
  resultsFreq2012_4 <- c(NA,NA,NA,NA,resultsFreq2012_4$freq_1)
  
  resultsFreq2012_5 <- tail(resultsFreq2012,-5) %>%
    mutate(freq_1 = n / sum(n)) 
  resultsFreq2012_5 <- c(NA,NA,NA,NA, NA, resultsFreq2012_5$freq_1)
  
  resultsFreq2012_6 <- tail(resultsFreq2012,-6) %>%
    mutate(freq_1 = n / sum(n)) 
  resultsFreq2012_6 <- c(NA,NA,NA,NA,NA,NA,resultsFreq2012_6$freq_1)
  
  resultsFreq2012_7 <- tail(resultsFreq2012,-7) %>%
    mutate(freq_1 = n / sum(n)) 
  resultsFreq2012_7 <- c(NA,NA,NA,NA,NA,NA,NA,resultsFreq2012_7$freq_1)
  
  
  resultsFreq2012 <- resultsFreq2012 %>%
    mutate(freq_1 = resultsFreq2012_1) %>%
    mutate(freq_2 = resultsFreq2012_2) %>%
    mutate(freq_3 = resultsFreq2012_3) %>%
    mutate(freq_4 = resultsFreq2012_4) %>%
    mutate(freq_5 = resultsFreq2012_5) %>%
    mutate(freq_6 = resultsFreq2012_6) %>% 
    mutate(freq_7 = resultsFreq2012_7) 
  
  resultsFreq2012
}

process2013DatSA <- function(resultsFreq2012, resultsFreq2013){
  
  results50 <- resultsFreq2013[resultsFreq2013$seroMax == 50,]$n * resultsFreq2012$freq
  results50 <- c(results50[1:7])
  results100 <- resultsFreq2013[resultsFreq2013$seroMax == 100,]$n * resultsFreq2012$freq_1
  results100 <- c(results100[2:7],NA)
  results200 <- resultsFreq2013[resultsFreq2013$seroMax == 200,]$n * resultsFreq2012$freq_2
  results200 <- c(results200[3:7],NA,NA)
  results400 <- resultsFreq2013[resultsFreq2013$seroMax == 400,]$n * resultsFreq2012$freq_3
  results400 <- c(results400[4:7],NA,NA,NA)
  results800 <- resultsFreq2013[resultsFreq2013$seroMax == 800,]$n * resultsFreq2012$freq_4
  results800 <- c(results800[5:7],NA,NA,NA,NA)
  results1600 <- resultsFreq2013[resultsFreq2013$seroMax == 1600,]$n * resultsFreq2012$freq_5
  results1600 <- c(results1600[6:7],NA,NA,NA,NA,NA)
  results3200 <- resultsFreq2013[resultsFreq2013$seroMax == 3200,]$n * resultsFreq2012$freq_6
  results3200 <- c(results3200[7],NA,NA,NA,NA,NA,NA)
  
  titreDropYear <- 0.66
  
  newDistribution <- matrix(nrow = 7, ncol = 8)
  colnames(newDistribution) <- c("50","100","200","400","800","1600","3200","years")
  newDistribution[,1] <- results50 
  newDistribution[,2] <- results100
  newDistribution[,3] <- results200
  newDistribution[,4] <- results400
  newDistribution[,5] <- results800
  newDistribution[,6] <- results1600
  newDistribution[,7] <- results3200
  newDistribution[,8] <- c(0,1*titreDropYear,2*titreDropYear,3*titreDropYear,4*titreDropYear,5*titreDropYear,6*titreDropYear)
  
  newDistribution <- tbl_df(newDistribution)
  
  newTallyDistribution <- newDistribution %>%
    replace(., is.na(.), 0) %>%
    mutate(yearsSum = rowSums(.[1:7])) %>%
    mutate(years = as.factor(years))
  
  totalYearSum <- sum(newTallyDistribution$yearsSum)
  
  newTallyDistribution <- newTallyDistribution %>%
    mutate(freq = yearsSum/totalYearSum)
  newTallyDistribution
}