## load the nimble library and set seed
load("linearGaussianSSM/simulatedDataEx1.RData")

# Load packages
library('nimble')
library(nimbleSMC)
library(nimMCMCSMCupdates)

#Type of particle filter to run
pfTypeRun = "bootstrap"

# load data
source("linearGaussianSSM/functionSimEstimation.R")


library(parallel)
library(doParallel)
thisCluster <- makeCluster(5)

# Fit the bootstrap PF
# We run the code in batches of 15 on a server
bootstrapEstimates <- parallel::parLapply(cl = thisCluster,
                                          X = 1:15,
                                          fun = runFunction,
                                          simData = simData,
                                          iNodePrev = c(49, 45, 20, 10, 5),
                                          nIterations = 30000,
                                          nBurnin = 20000,
                                          nChains = 2,
                                          nThin= 1,
                                          nyears = 50,
                                          numParticles = 1000,
                                          pfTypeRun = "bootstrap")



#save results
save(bootstrapEstimates, file = "linearGaussianSSM/bootstrapPF/estimates1.RData")
stopCluster(thisCluster)

