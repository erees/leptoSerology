
mcmcModelRandomSampler <- function(numberOfSamples, mcmcMatrix, ageVector){
  for (i in 1:numberOfSamples){
    # Specify random number to sample from mcmc chain for parameter estimates
    randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
    lambdaSample <- mcmcMatrix[randomNumber,2]
    deltaSample <- mcmcMatrix[randomNumber,1]
    # Run the model with new esitmates
    newRow <- (lambdaSample / (lambdaSample+deltaSample)) *(1 - exp(-ageVector*(lambdaSample+deltaSample)))
    # Add this to a new matrix, where each row is a sample, and each column is an age
    outDf[i,] <- newRow
  }
  outDf
}

mcmcRandomSamplerCat <- function(numberOfSamples, mcmcMatrix, ageVector, ageTotals){
  
  outDf <- matrix(NA,nrow=numberOfSamples, ncol = length(ageVector))
  
  
  for (i in 1:numberOfSamples){
    randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
    
    lambdaSample <- mcmcMatrix[randomNumber,1]
    
    newRow <- 1 - exp(-ageVector*(lambdaSample))
    updateRow <- c()
    for(j in 1:length(ageTotals)){
      randomlySampleBinomialDis <- rbinom(1,size = ageTotals[j],prob = newRow[j])
      if(randomlySampleBinomialDis > 0){
        result <- randomlySampleBinomialDis / ageTotals[j]
      } else {
        result <- 0
      }
      updateRow[j] <- result
    }
    rbinom(1,20,0.03)
    outDf[i,] <- updateRow
  }
  outDf
}


mcmcRandomSampler <- function(numberOfSamples, mcmcMatrix, ageVector, ageTotals){
  
  outDf <- matrix(NA,nrow=numberOfSamples, ncol = length(ageVector))
  
  
  for (i in 1:numberOfSamples){
    randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
    
    lambdaSample <- mcmcMatrix[randomNumber,2]
    deltaSample <- mcmcMatrix[randomNumber,1]
    
    newRow <- (lambdaSample / (lambdaSample+deltaSample)) *(1 - exp(-ageVector*(lambdaSample+deltaSample)))
    updateRow <- c()
    for(j in 1:length(ageTotals)){
      randomlySampleBinomialDis <- rbinom(1,size = ageTotals[j],prob = newRow[j])
      if(randomlySampleBinomialDis > 0){
        result <- randomlySampleBinomialDis / ageTotals[j]
      } else {
        result <- 0
      }
      updateRow[j] <- result
    }
    
    outDf[i,] <- updateRow
  }
  outDf
}

ageQuantiles <- function(mcmcDF){
  quantileMatrix <- matrix(NA,nrow=ncol(mcmcDF), ncol = 3)
  for(i in 1:ncol(mcmcDF)){
    quantiles <- mcmcDF[,i] %>% quantile(probs=c(.5,.025,.975))
    quantileMatrix[i,] <- quantiles
  }
  quantileMatrix
}
