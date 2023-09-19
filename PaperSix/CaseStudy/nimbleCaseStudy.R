#Load the required data
#load("CaseStudy/data_idm.RData")
#source("cross_valid_function.R")
#nimbleOptions(MCMCusePredictiveDependenciesInCalculations = FALSE)
#nimbleOptions(determinePredictiveNodesInModel = TRUE)
#nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)

# Packages required to load this script
require(coda)
require(nimble)
require(MCMCglmm)
library(parallel)
library(pbapply)
require(ggmcmc)
library(doParallel)
library(nimble)

op <- pboptions(type="timer")

# Function to estimate sum
#takes away infinite numbers
mysum <- function(x){
  ret <- sum(x[!is.infinite(x)], na.rm = TRUE)
  return(ret)
}

formatMatrix <- function(x, nSite, nVisit){
  newX <- c(x)
  y <- matrix(newX, nrow = nSite, ncol = nVisit)
  ret <- rowMeans(y)
  return(ret)
}

#compile the R function
nimble_sum <- nimbleRcall(
  prototype = function(
    x=double(1)
  ) {},
  returnType = double(0), # outcome is a number
  Rfun = 'mysum'
)

nimbleFormatMatrix <- nimbleRcall(
  prototype = function(
    x=double(1),
    nSite = double(0),
    nVisit = double(0)
  ) {},
  returnType = double(1), # outcome is a number
  Rfun = 'formatMatrix'
)

# Function to analyse the data and returns output for plots and summary
# Takes the dataformated from the dataFormat.R script
# method: whether IDM, Species only (spe) or group count (IG) only
# crossvalidate is a T/F indicating whether crossvalidation should be performed
# model refers to the IDM framework used: shared, covariate, correlation
# covariancePrior: "LV" for the multiplicative gamma shrinkage process and "full"
# for the covariance prior from the inverse Wishart distribution

run_nimble_model <- function(simulations_all,
                             method,
                             crossvalidate = FALSE,
                             model,
                             covariance_prior = "LV"){

  if(!method %in% c("covariate", "shared")) stop("Model can only be covariate or shared")
  if(!model %in% c("IDM", "Spe", "IG")) stop("model must be IDM, IG or Spe")
  if(!covariance_prior %in% c("LV", "full")) stop("Covariance prior must either be 'full' or 'LV'")
  if(!is.logical(crossvalidate)) stop("crossvalidate must be True or False")

  #functions for estimation in NIMBLE Code
  mysum <- function(x){
    ret <- sum(x[!is.infinite(x)], na.rm = TRUE)
    return(ret)
  }

  formatMatrix <- function(x, nSite, nVisit){
    newX <- c(x)
    y <- matrix(newX, nrow = nSite, ncol = nVisit)
    ret <- rowMeans(y)
    return(ret)
  }

  #compile the R function
  nimble_sum <- nimbleRcall(
    prototype = function(
    x=double(1)
    ) {},
    returnType = double(0), # outcome is a number
    Rfun = 'mysum'
  )

  nimbleFormatMatrix <- nimbleRcall(
    prototype = function(
    x=double(1),
    nSite = double(0),
    nVisit = double(0)
    ) {},
    returnType = double(1), # outcome is a number
    Rfun = 'formatMatrix'
  )


   start_time <- Sys.time()

  #write nimble code
  code <- nimbleCode({

    #############################################################
    #                 PRIOR DISTRIBUTIONS                       #
    #############################################################

    #specific specific intercept for species occupancy model
    #same for the group count model using the shared model
    for(spe.tag in 1:n.species){
      betaSpecies[spe.tag]~ dnorm(muSpecies, sd=sigmaSpecies)
      betaLatitude[spe.tag] ~ dnorm(muLatitude, sd = sigmaLatitude)
    }

    # Intercept to be used for the covariate structure for the group count data
    betaLambda ~ dnorm(mulambda, sd = sigmaLambda)

    #overdispersion parameter
    rho ~ dunif(0.01, 50)

    #species interaction mean
    # for(spe.tag in 1:n.species){
    #   mu.eta[spe.tag] <- 0
    # }

    #parameter priors
   mulambda ~ dnorm(0, sd = 10)
    muSpecies ~ dnorm(0, sd = 10) #community intercept
    muLatitude ~ dnorm(0, sd = 10) #community slope of covariate
    muAlphaSpecies <- 0 #detection prob
    muAlphaSites <- 0 #detection prob
    muVisits <- 0 #detection prob
    sigmaSpecies ~ T(dnorm(0, 0.1), 0.1, 10) #community intercept standard deviation
    sigmaLatitude ~ T(dnorm(0, 0.1), 0.1, 10) #community slope standard deviation
    sigmaLambda ~ T(dnorm(0, 0.1), 0.1, 10) # standard deviation of intercept of GCM under the covariate model
    sigmaAlphaSpecies ~ T(dnorm(0, 0.1), 0.1, 10) #detection prob species random effect standard deviation
    sigmaAlphaSites~ T(dnorm(0, 0.1), 0.1, 10) #detection prob sites random effect standard deviation
    #sigmaVisits~ T(dnorm(0, 0.1), 0.1, 5) #detection prob visits random effect standard deviation
    sigmaSites ~ T(dnorm(0, 0.1), 0.1, 10)
    # detection process species effect
    for(spe.tag in 1:n.species){
      alphaSpecies[spe.tag]~ dnorm(muAlphaSpecies, sd = sigmaAlphaSpecies)
    }

    # detection process sites effect
    for(site.tag in 1:n.sites){
      alphaSites[site.tag]~ dnorm(muAlphaSites, sd = sigmaAlphaSites)
    }

    # detection process visit effect
    # for(visit.tag in 1:n.visit){
    #   betaVisits[visit.tag] ~ dnorm(muVisits, sd = sigmaVisits)
    # }

    #Detection probability
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
       # for(visit.tag in 1:n.visit){
          logit(p.tag[site.tag, spe.tag]) <- alphaSites[site.tag] + alphaSpecies[spe.tag] #+  betaVisits[visit.tag]
       #}
      }
    }


    # Site effect
    for(site.tag in 1:n.sites){
      eta.lam[site.tag] ~ dnorm(0, sd = sigmaSites)
    }
    #################################
    # covariance prior for interraction effect
    ##############################
    # if(covariance_prior == "full"){
    #   omega[1:n.species,1:n.species] ~ dwish(R[1:n.species,1:n.species], df)
    #   Cov[1:n.species,1:n.species] <- inverse(omega[1:n.species,1:n.species])
    # }

    # if(covariance_prior == "LV"){
    #   delta[1] ~ dgamma(a1,1)
    #   for(factor in 2:NLatent) {
    #     delta[factor] ~ dgamma(a2,1)
    #   }
    #
    #   for(l in 1:NLatent){
    #     tauDelta[l] <- prod(delta[1:l])
    #   }
    #
    #   for(spe.tag in 1:n.species){
    #     sig[spe.tag] ~ dgamma(a.sigma, b.sigma)
    #   }
    #   Sigma[1:n.species, 1:n.species] <- diag(1/sig[1:n.species])
    #
    #   for (spe.tag in 1:n.species) {
    #     for (l in 1:NLatent) {
    #       phi[spe.tag,l] ~ dgamma(nu/2,nu/2)
    #       tauFS[spe.tag, l] <- 1/(phi[spe.tag,l]*tauDelta[l])
    #       lamLatent[spe.tag, l] ~ dnorm(mean = 0, sd = sqrt(tauFS[spe.tag,l]))
    #     }
    #   }
    #   Cov[1:n.species, 1:n.species] <- lamLatent[1:n.species, 1:NLatent] %*% t(lamLatent[1:n.species, 1:NLatent]) + Sigma[1:n.species, 1:n.species]
    #   omega[1:n.species,1:n.species] <- inverse(Cov[1:n.species, 1:n.species])
    # }


    # for(site.tag in 1:n.sites){
    #   # species interraction effect
    #   eta.lam[site.tag, 1:n.species] ~ dmnorm(mean = mu.eta[1:n.species],
    #                                           prec = omega[1:n.species, 1:n.species])
    # }

    #correlation matrix
    # for (spe.tag in 1:n.species) {
    #   for (ss in 1:n.species) {
    #     CorrIn[spe.tag, ss] <- Cov[spe.tag, ss]/sqrt(Cov[spe.tag, spe.tag] * Cov[ss, ss])
    #   }
    # }

    ####################################
    #Link between the abundance and occupancy -
    # Defining linear predictors
    #####################################
    if(method == "covariate"){
      for(site.tag in 1:n.sites){ #loop over sites
        for(spe.tag in 1:n.species){#loop over species
          #for(visit.tag in 1:n.visit){ #loop over visits
          muAll[site.tag,spe.tag] <- betaSpecies[spe.tag] + eta.lam[site.tag] + betaLatitude[spe.tag] * latitude[site.tag]  + betaLambda
          muPsi[site.tag,spe.tag] <- betaSpecies[spe.tag] + eta.lam[site.tag] + betaLatitude[spe.tag] * latitude[site.tag]
          muLambda[site.tag,spe.tag] <- eta.lam[site.tag] + betaLatitude[spe.tag] * latitude[site.tag] + betaLambda
          log(lambda[site.tag, spe.tag]) <- muLambda[site.tag, spe.tag]
          cloglog(psi[site.tag, spe.tag]) <- muPsi[site.tag, spe.tag]
         # }
        }
      }
    }

    if(method == "shared"){
      for(site.tag in 1:n.sites){ #loop over sites
        for(spe.tag in 1:n.species){#loop over species
         # for(visit.tag in 1:n.visit){
          muAll[site.tag,spe.tag] <- betaSpecies[spe.tag] + eta.lam[site.tag] + betaLatitude[spe.tag] * latitude[site.tag]
          muPsi[site.tag,spe.tag] <- betaSpecies[spe.tag] + eta.lam[site.tag] + betaLatitude[spe.tag] * latitude[site.tag]
          muLambda[site.tag,spe.tag] <- eta.lam[site.tag] + betaLatitude[spe.tag] * latitude[site.tag] + betaSpecies[spe.tag]
          log(lambda[site.tag, spe.tag]) <- muLambda[site.tag,spe.tag]
          cloglog(psi[site.tag, spe.tag]) <- muPsi[site.tag,spe.tag]
         # }
        }
      }
    }


    ###############################################
    #   JOINT LIKELIHOOD FOR INTEGRATED DISTRITION MODELS
    ################################################

    if(model == "IDM"){
      #Model for Group count data
      for(site.tag in 1:n.sites){ # loop over sites
        lambda.g[site.tag] <- sum(lambda[site.tag,1:n.species])
        lambda.p[site.tag] <- rho /(rho + lambda.g[site.tag])
        for(visit.tag in 1:n.visit){ #loop over visits
           Y[site.tag, visit.tag] ~ dnegbin(lambda.p[site.tag], rho)
        }
      }

      #   #Model for species occupancy data
       for(site.tag in 1:n.sites){ # loop over sites
         for(spe.tag in 1:n.species){ #loop over species
           #for(visit.tag in 1:n.visit){ #loop over visits
             z[site.tag, spe.tag] ~ dbern(psi[site.tag, spe.tag])
          # }
         }
       }

      for(site.tag in 1:n.sites){ # loop over sites
        for(spe.tag in 1:n.species){ #loop over species
          for(visit.tag in 1:n.visit){ #loop over visits
           # z[site.tag, spe.tag, visit.tag] ~ dbern(psi[site.tag, spe.tag, visit.tag])
            X[site.tag, spe.tag, visit.tag] ~ dbin(z[site.tag, spe.tag]*p.tag[site.tag,spe.tag], n.replicates)
          }
        }
      }

      # Relative proportion for diversity indices
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          #for(visit.tag in 1:n.visit){
          pps[site.tag, spe.tag] <- exp(muAll[site.tag, spe.tag])
          pis[site.tag, spe.tag] <- pps[site.tag, spe.tag]/sum(pps[site.tag, 1:n.species])
         # }
        }
      }
    }

    ###############################################
    #   LIKELIHOOD FOR GROUP COUNT  MODELS
    ################################################
     if(model == "IG"){

       for(site.tag in 1:n.sites){ # loop over sites
         for(spe.tag in 1:n.species){ #loop over species
           #for(visit.tag in 1:n.visit){ #loop over visits
           z[site.tag, spe.tag] ~ dbern(psi[site.tag,spe.tag])
           #   X[site.tag,spe.tag,visit.tag] ~ dbin(z[site.tag, spe.tag, visit.tag]*p.tag[site.tag,spe.tag], n.replicates)
           #}
         }
         #}
       }


       for(site.tag in 1:n.sites){ # loop over sites

         lambda.p[site.tag] <- rho /(rho + lambda.g[site.tag])
         lambda.g[site.tag] <- sum(lambda[site.tag,1:n.species])
         for(visit.tag in 1:n.visit){ #loop over visits
           Y[site.tag, visit.tag] ~ dnegbin(lambda.p[site.tag], rho)
         }
       }
    #   # Relative proportion for diversity indices
      for(site.tag in 1:n.sites){
         for(spe.tag in 1:n.species){
           pps[site.tag, spe.tag] <- exp(muLambda[site.tag, spe.tag])
          pis[site.tag, spe.tag] <- pps[site.tag, spe.tag]/sum(pps[site.tag, 1:n.species])
         }
       }
    }


    ###############################################
    #   LIKELIHOOD FOR SPECIES OCCUPANCY MODELS
    ################################################
    if(model == "Spe"){
      for(site.tag in 1:n.sites){ # loop over sites
        for(spe.tag in 1:n.species){ #loop over species
          #for(visit.tag in 1:n.visit){ #loop over visits
            z[site.tag, spe.tag] ~ dbern(psi[site.tag,spe.tag])
         #   X[site.tag,spe.tag,visit.tag] ~ dbin(z[site.tag, spe.tag, visit.tag]*p.tag[site.tag,spe.tag], n.replicates)
          #}
        }
        #}
      }

      for(site.tag in 1:n.sites){ # loop over sites
        for(spe.tag in 1:n.species){ #loop over species
          for(visit.tag in 1:n.visit){ #loop over visits
            X[site.tag,spe.tag, visit.tag] ~ dbin(z[site.tag, spe.tag]*p.tag[site.tag,spe.tag], n.replicates)
          }
        }
      }

      # Relative proportion for diversity indices
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
         # for(visit.tag in 1:n.visit){ #loop over visits
          pps[site.tag, spe.tag] <- exp(muPsi[site.tag, spe.tag])
          pis[site.tag, spe.tag] <- pps[site.tag, spe.tag]/sum(pps[site.tag, 1:n.species])
         # }
          }
      }
    }

    #     #############################################################
    #     #               shannon Index                              #
    #############################################################

    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        #for(visit.tag in 1:n.visit){
        shans[site.tag, spe.tag] <- log(pis[site.tag, spe.tag])*(pis[site.tag, spe.tag])
      #}
      }
    }

    for(site.tag in 1:n.sites){
      #for(visit.tag in 1:n.visit){
      shan[site.tag] <- - nimble_sum(shans[site.tag, 1:n.species])
      #}
    }

    #for(site.tag in 1:n.sites){
      #averageShan[site.tag] <- sum(shan[site.tag, 1:n.visit])/n.visit
    #}
    #averageShan[1:(n.sites/8)] <- nimbleFormatMatrix(shan[1:n.sites], (n.sites/8), 8 )
  })

  ###########################
  # DETAILS FOR THE NIMBLE MCMC
  ##########################

  # Extracting dimensions of data parameters
  data <- simulations_all

  dim_data <- dim(data[[1]])
  #data_dim <- dim_data[2]

  #constants for the model
  const <- list(n.sites = dim_data[1],
                n.species= (dim_data[2]),
               # n.years = dim_data[4],
                n.replicates = 5,
               visit= data$visit,
               n.visit = length(unique(data$visit)),
                mu.eta = rep(0, dim_data[2]),
                nu = 3,
                NLatent = ceiling((dim_data[2])/5),
                a1 = ceiling((dim_data[2])/2),
                a2= ceiling((dim_data[2])/2) + 5,
                a.sigma = 1,
                b.sigma = 1
  )

  #R and df for the inverse Wishart prior
  Rmat <- diag(const$n.species) # In the paper, this is initial covariance matrix
  df <- const$n.species + 1

  latitude <- simulations_all[[4]]
  latitude[is.na(latitude)] <- 100000
  latitude <- as.numeric(scale(latitude, center = TRUE, scale = TRUE))

  pa_data <- function(data){
    zst <- apply(data, c(1,2), function(x) max(x, na.rm = TRUE))
    #zst <- data
    zst[zst!= 0] <- 1
    zst[zst == 0] <- NA
    #zst[zst!= 0] <- 1
    return(zst)
  }
  zst = pa_data(simulations_all[[1]])

  #Data for the model
  idm_data <- list(Y = round(simulations_all[[2]]/data$Nsurveys),
                   X = simulations_all[[1]],
                   latitude = rep(latitude),
                   #Nsurveys = log(data$Nsurveys),
                   R = Rmat,
                   df = df,
                   z = zst)



  #For LV approach
  delta = rgamma(const$NLatent, 2,1)
  phi =matrix(rgamma(const$n.species*const$NLatent,1.5,1.5), ncol=const$NLatent, nrow=const$n.species)
  tau= cumprod(delta)
  lamLatent = matrix(NA, ncol=const$NLatent, nrow=const$n.species)
  for(spe.tag in 1:const$n.species){
    for(i in 1:const$NLatent){
      lamLatent[spe.tag,i] <- rnorm(1,0,sd=sqrt(1/(phi[spe.tag,i]*tau[i])))
    }
  }

  # Initial values for the model
  idm_inits <- function(){list(betaSpecies= rnorm(const$n.species,0,0.1),
                               betaLatitude= rnorm(const$n.species,0,0.1),#covariate effect
                               betaLambda = rnorm(1,0,0.1),
                               alphaSpecies=rnorm(const$n.species,0,0.1), #detection probability
                               alphaSites=rnorm(const$n.sites,0,0.1),
                               betaVisits=rnorm(const$n.visit,0,0.1),
                               sigmaSpecies = 1,
                               sigmaLatitude = 1,
                               sigmaLambda = 1,
                               sigmaAlphaSpecies = 1,
                               sigmaAlphaSites = 1,
                               sigmaVisits = 1,
                               sigmaSites = 1,
                               muSpecies = 0,
                               muLatitude = 0,
                               muVisits = 0,
                               #z= zst,
                               #delta=delta,
                               #phi=phi,
                               rho = 0.8,
                               mulambda = 0,
                               #lamLatent= lamLatent ,
                               #omega=diag(1, const$n.species),
                               #sig = rgamma(const$n.species, 2,1),
                               eta.lam = rep(0, const$n.sites)
  )
  }

  initsList <- idm_inits()

  #############################################################
  #                 Compile Model                             #
  #############################################################

  mwtc <- nimbleModel(code,
                      name = "mwtc",
                      data = idm_data,
                      constants = const,
                      inits = initsList)

  Cmwtc <- compileNimble(mwtc)

  mcmcconf <- configureMCMC(Cmwtc,
                            #print=FALSE,
                            monitors = c("muSpecies",
                                         "muLatitude",
                                         #"muVisits",
                                         #"sigmaVisits",
                                         "sigmaSpecies",
                                         "sigmaLatitude",
                                         "sigmaAlphaSpecies",
                                         "sigmaAlphaSites",
                                         "sigmaLambda",
                                         "sigmaSites",
                                         "rho",
                                         #"Cov",
                                         "eta.lam",
                                         "shan",
                                         "betaSpecies",
                                         "betaLatitude",
                                         "betaLambda",
                                         #"betaVisits",
                                         "alphaSpecies",
                                         "alphaSites",
                                         "mulambda",
                                         "z"
                                         #"lamLatent",
                                         #"phi",
                                         #"sig",
                                         #"delta",
                                         #"betaLambda",
                                        # "z"
                                         ),
                            enableWAIC = TRUE)
# mcmcconf$removeSamplers(c("betaSpecies",
#                           "betaLatitude",
#                           "alphaSpecies",
#                           "alphaSites"), print = FALSE)
  #mcmcconf$addSampler(target = "delta", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "sig", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "phi", type = "RW_block", print = TRUE)
 # mcmcconf$addSampler(target = "lamLatent", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "betaSpecies", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "betaLatitude", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "alphaSpecies", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "alphaSites", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "betaVisits", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "betaVisits", type = "AF_slice", print = TRUE)
#build MCMC and compile
  Rmcmc <- buildMCMC(mcmcconf)


  cmcmc <- compileNimble(Rmcmc,
                         project = Cmwtc)

  #run MCMC
  start_time <- Sys.time()
  estimate <- runMCMC(cmcmc,
                      niter = 50000,
                      nburnin = 30000,#100000
                      inits = initsList,
                      nchains = 3,
                      thin = 4, #one
                      summary=TRUE,
                      samples=TRUE,
                      setSeed = TRUE,
                      samplesAsCodaMCMC=TRUE,
                      WAIC = TRUE)

  end.time <- Sys.time()
  time_taken <- as.numeric(round(end.time-start_time,2))

  estimation <- estimate$summary$all.chains

  ret <- NULL

  return(list(estimation,time_taken, estimate, ret))

}




