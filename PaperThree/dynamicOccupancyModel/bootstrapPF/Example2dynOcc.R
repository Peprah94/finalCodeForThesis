# This script is used to run the dynamic occupancy model described
# the second simulation study using the boostrap PF

# Load packages
library(nimble)
library(nimbleSMC)
library(nimMCMCSMCupdates)

# set the configurations for fitting the model
nyears = 55
nsites = 300
nvisits = 6
iNodePrev = 45
pfTypeRun = "bootstrap"

# Set the configurations for MCMC
nIterations = 100000
nBurnin = 50000
nChains = 2
nThin = 10
numParticles = 50

# Load simulated data and results from reduced model from auxiliary PF folder
load("simDataDynamicOccupancy.RData")

# NIMBLE code
dynOccupancyCode <- nimbleCode({

  # Prior distributions of hyperparameters
  alphaPSig ~ T(dnorm(0, 0.1), 0.001, )
  betaPSig ~ T(dnorm(0, 0.1), 0.001, )
  alphaPsiSig ~ T(dnorm(0, 0.1), 0.001, )
  betaPsiSig ~ T(dnorm(0, 0.1), 0.001, )


  alphaPhi ~ dnorm(0, sd = 10)
  betaPhi ~ dnorm(0, sd = 10)

  # Prior distributions
    alphaP ~ dnorm(mean = 0, sd = alphaPSig)
    betaP ~ dnorm(mean = 0, sd =  betaPSig)
    alphaPsi ~ dnorm(mean = 0, sd =  alphaPsiSig)
    betaPsi ~ dnorm(mean = 0, sd =  betaPsiSig)

  # Detection Probability
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        logit(detectProb[site.tag, visit.tag, year.tag]) <- alphaP + betaP* windSpeed[site.tag, visit.tag, year.tag]
      }
    }
  }

  # Initial occupancy probability psi
  for(site.tag in 1:nsites){
    logit(initOccuProb[site.tag]) <- alphaPhi + betaPhi*elevation[site.tag]
  }

  # Persistence and colonisation probability
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      logit(persistenceProb[site.tag,  year.tag]) <- alphaPsi + betaPsi* springPrecipitation[site.tag, year.tag]
      # logit(colonisationProb[site.tag,  year.tag]) <- alphaGamma[year.tag] + betaGamma[year.tag]* sizeOfBeak[site.tag, year.tag]
    }
  }

    colonisationProb ~ dunif(0.01, 0.5)


  # Initial Presence/Absence
  for(site.tag in 1:nsites){
    z[site.tag, 1] ~ dbin(prob = initOccuProb[site.tag], 1)
  }

  # True presence/absence
  for(year.tag in 2:nyears){
    for(site.tag in 1:nsites){
      z[site.tag, year.tag] ~ dbin(prob = (z[site.tag, (year.tag -1)] * persistenceProb[site.tag,  year.tag] + (1 - z[site.tag, (year.tag -1)])*colonisationProb), 1)
    }
  }

  #observations
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        y[site.tag, visit.tag, year.tag] ~ dbin(prob = z[site.tag, year.tag] * detectProb[site.tag, visit.tag, year.tag], 1)
      }
    }
  }

  # Derived quantities
  for(year.tag in 1:nyears){
    psi.fs[year.tag] <- sum(z[1:nsites, year.tag])/nsites
  }

})


# ## define data, constants, and initial values
z <- apply(simData$y, c(1, 3), function(x){
  ifelse(sum(x) > 0, 1, NA)
})


#####################
#   Reduced Model
########################

data <- list(
  y = simData$y[,,-c((iNodePrev +1):nyears)],
  windSpeed = simData$covariates$windSpeed[,,-c((iNodePrev +1):nyears)],
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation[,-c((iNodePrev +1):nyears)],
  sizeOfBeak = simData$covariates$sizeOfBeak,
  z = z[,-c((iNodePrev +1):nyears)]
)

constants <- list(
  nyears = iNodePrev,
  nsites = nsites,
  nvisits = nvisits
)
inits <- list(
  alphaPSig = 1,
  betaPSig  = 1,
  alphaPsiSig = 1,
  betaPsiSig  = 1,
  alphaGammaSig = 1,
  betaGammaSig = 1,
  alphaPhi  = 0,
  betaPhi  = 0,
  alphaP =  rnorm(1, mean = 0, sd = 1),
  betaP =  rnorm(1, mean = 0, sd = 1),
  alphaPsi=  rnorm(1, mean = 0, sd = 1),
  betaPsi= rnorm(1, mean = 0, sd = 1),
  colonisationProb = 0.01,
  alphaPhi=  rnorm(1),
  betaPhi= rnorm(1)
)

# NIMBLE model for reduced model
newModelReduced <- nimbleModel(dynOccupancyCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)

# Copy results from auxiliary PF folder here.
load("example4ReducedBootstrapTrue.RData")

################
# Updated Model
################
data <- list(
  y = simData$y,
  windSpeed = simData$covariates$windSpeed,
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation,
  sizeOfBeak = simData$covariates$sizeOfBeak,
  z = z
)

constants <- list(
  nyears = nyears,
  nsites = nsites,
  nvisits = nvisits
)

inits <- list(
  alphaPSig = 1,
  betaPSig  = 1,
  alphaPsiSig = 1,
  betaPsiSig  = 1,
  alphaGammaSig = 1,
  betaGammaSig = 1,
  alphaPhi  = 0,
  betaPhi  = 0,
  alphaP =  rnorm(1, mean = 0, sd = 1),
  betaP =  rnorm(1, mean = 0, sd = 1),
  alphaPsi=  rnorm(1, mean = 0, sd = 1),
  betaPsi= rnorm(1, mean = 0, sd = 1),
  colonisationProb = 0.01,
  alphaPhi=  rnorm(1),
  betaPhi= rnorm(1)
)

# NIMBLE model for updated model
newModelUpdated <- nimbleModel(dynOccupancyCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)

# Fitting the model
example2UpdatedModelTrue <- spartaNimUpdates(model = newModelUpdated, 
                                             reducedModel = newModelReduced,
                                             latent = "z",
                                             nParFiltRun = numParticles,
                                             pfType = pfTypeRun,
                                             propCov = c(1,1,1,1,1,1,0.01)*diag(7),
                                             mcmcScale = 1,
                                             extraVars = NULL,
                                               MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
                                                                                 'alphaPsiSig', 'betaPsiSig',
                                                                                 'alphaPhi', 'betaPhi',
                                                                                 'alphaP', 'betaP',
                                                                                 'alphaPsi', 'betaPsi',
                                                                                 "colonisationProb"),
                                                                      additionalPars = c("z", "psi.fs"),
                                                                      n.iter = (nIterations - nBurnin)/nThin,
                                                                      n.chains = nChains,
                                                                      n.burnin = 10,
                                                                      n.thin = 1
                                             ),  #saved loglikelihoods from reduced model
                                             postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
                                             pfControl = list(saveAll = TRUE,
                                                              timeIndex = 2,
                                                              smoothing = TRUE,
                                                              mcmc = TRUE,
                                                              M = nyears - iNodePrev,
                                                              iNodePrev = iNodePrev)
)

save(example2UpdatedModelTrue, file= "example4UpdatedBootstrapTrue1.RData")


#############
# Baseline Model
###############
z <- apply(simData$y, c(1, 3), function(x){
  ifelse(sum(x) > 0, 1, NA)
})

data <- list(
  y = simData$y,
  windSpeed = simData$covariates$windSpeed,
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation,
  sizeOfBeak = simData$covariates$sizeOfBeak,
  z = z
)

constants <- list(
  nyears = nyears,
  nsites = nsites,
  nvisits = nvisits
)
inits <- list(
  alphaPSig = 1,
  betaPSig  = 1,
  alphaPsiSig = 1,
  betaPsiSig  = 1,
  alphaGammaSig = 1,
  betaGammaSig = 1,
  alphaPhi  = 0,
  betaPhi  = 0,
  alphaP =  rnorm(1, mean = 0, sd = 1),
  betaP =  rnorm(1, mean = 0, sd = 1),
  alphaPsi=  rnorm(1, mean = 0, sd = 1),
  betaPsi= rnorm(1, mean = 0, sd = 1),
  alphaPhi=  rnorm(1),
  betaPhi= rnorm(1)
)
#
#
# ## build the model
dynOccModel <- nimbleModel(dynOccupancyCode,
                           data = data,
                           constants = constants,
                           inits = inits,
                           check = FALSE)


# Fitting Baseline model
newModel <- dynOccModel$newModel(replicate = TRUE)
