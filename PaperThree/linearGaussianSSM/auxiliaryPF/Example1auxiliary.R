
## load the nimble library and set seed
load("linearGaussianSSM/simulatedDataEx1.RData")

# Load packages
library(nimble)
library(nimbleSMC)
library(nimMCMCSMCupdates)

#Particle filter to run
pfTypeRun = "auxiliary"

# load function to fit the models with MCMC
source("linearGaussianSSM/functionSimEstimation.R")


library(parallel)
library(doParallel)
thisCluster <- makeCluster(5)

# Fit an auxiliary PF
#We run code in batches of 15
auxiliaryEstimates <- parallel::parLapply(cl = thisCluster,
                                          X = 31:45,
                                          fun = runFunction,
                                          simData = simData,
                                          iNodePrev = c(49, 45, 20, 10, 5),
                                          nIterations = 30000,
                                          nBurnin = 20000,
                                          nChains = 2,
                                          nThin= 1,
                                          nyears = 50,
                                          numParticles = 1000,
                                          pfTypeRun = "auxiliary")


#save results
save(auxiliaryEstimates, file = "linearGaussianSSM/auxiliaryPF/estimates3.RData")
stopCluster(thisCluster)

