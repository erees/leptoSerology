extractLogLik <- function(FOIfit){
  
  estimated_parameters <- FOIfit$model$estimated_parameters
  chains <- rstan::extract(FOIfit$fit)
  FOIs <- chains$lambda
  bg <- chains$bg
  M <- nrow(FOIs)
  N <- FOIfit$data$N
  A <- FOIfit$data$A
  Ncategory <- FOIfit$data$Ncategory
  NAgeGroups <- FOIfit$data$NAgeGroups
  LogLikelihoods <- matrix(0, nrow = M, ncol = N)
  Y <- FOIfit$data$Y
  age <- FOIfit$data$age
  category <- FOIfit$data$categoryindex
  for (i in seq(1, M)) {
    lk = chains$Like[i, ]
    for (j in seq(1, N)) {
      if (Y[j] == FALSE) {
        L = log(1 - lk[j])
      }
      else {
        L = log(lk[j])
      }
      LogLikelihoods[i, j] <- L
    }
  }
  return(LogLikelihoods)
}