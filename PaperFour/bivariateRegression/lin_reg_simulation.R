
#Packages needed to run the model
library(nimble)
library(INLA)
library(mvtnorm)
library(MASS)
library(parallel)
library(coda)
library(ggmcmc)

# function for generating samples
set.seed(0)
sample.linreg <- function(){
  n = 100 # Number of samples
  x1 = runif(n) #covariate 1
  x2 = runif(n) #covariate 2
  err = rnorm(n) # error
  y = 2 + 3*x1 -3*x2 + err # Response
  return(list(y = y,x = as.matrix(cbind(x1, x2))))
}

#Data
bivariateSims <- sample.linreg() #Data from the sampling fnx
save(bivariateSims, file="data_for_simulations.RData")















