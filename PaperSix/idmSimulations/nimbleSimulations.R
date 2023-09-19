#Runs the analyses of the simulated data
#takes input the simulated data, methods ("IDM", "Spe", "IG"), and
#covariance_prior ("full" for sigma from invWishart), and "LV" for multiplicative
#shrinkage prior (a form of latent variable approach for sigma)
#and outputs the summary of the alpha's, beta's, z's, lambda's and correlation matrix
#as well as the values of alpha's and beta's that have converged


#Packages required to run the package
library(pbapply)
library(nimble)
library(MCMCglmm)
library(coda)
library(MCMCpack)
library(boot)
library(MASS)
library(parallel)
require(ggmcmc)


# I defined the functions here for the sake of parallelisation
incidence <- function(lambdas, prob, K=4){ #K is the number of visits
  incidence_mat <- lambdas*(1-(1-prob)^K)
  return(incidence_mat)
}

nimble_incidence <- nimbleRcall(
  prototype = function(
    lambda=double(2),
    prob = double(2),
    K = double(0)
    # beta is a vector
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'incidence'
)

richness <- function(z){
  rowSums(z)
}

nimble_richness <- nimbleRcall(
  prototype = function(
    z= double(2)
  ) {},
  returnType = double(1), # outcome is a vector
  Rfun = 'richness'
)

hill_index <- function(q, incidence){
  pis <- proportions(incidence, margin = 1)
  if(q != 1){
    hill <- (rowSums(pis^q, na.rm = TRUE))^(1/(1-q))
  }else{
    hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))
  }
  return(hill)
}

nimble_hill_index <- nimbleRcall(
  prototype = function(
    q=double(0),
    incidence = double(2)
  ) {},
  returnType = double(1), # outcome is a vector
  Rfun = 'hill_index'
)

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

run_nimble_model <- function(simulations_all, method, covariance_prior = "LV", model){
  library(pbapply)
  library(nimble)
  library(MCMCglmm)
  library(coda)
  library(MCMCpack)
  library(boot)
  library(MASS)
  library(parallel)
  require(ggmcmc)
  #checks
  if(!method %in% c("covariate", "shared")) stop("Model can only be covariate or shared")
  if(!model %in% c("IDM", "Spe", "IG")) stop("model must be IDM, IG or Spe")
  if(!covariance_prior %in% c("LV", "full")) stop("Covariance prior must either be 'full' or 'LV'")
  #if(!is.logical(crossvalidate)) stop("crossvalidate must be True or False")

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
    for(spe.tag in 1:n.species){
      mu.eta[spe.tag] <- 0
    }

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
    sigmaAlphaSites ~ T(dnorm(0, 0.1), 0.1, 10) #detection prob sites random effect standard deviation
    sigmaSites ~ T(dnorm(0, 0.1), 0.1, 10) #detection prob visits random effect standard deviation

    # detection process species effect
    for(spe.tag in 1:n.species){
      alphaSpecies[spe.tag]~ dnorm(muAlphaSpecies, sd = sigmaAlphaSpecies)
    }

    # detection process sites effect
    for(site.tag in 1:n.sites){
      alphaSites[site.tag]~ dnorm(muAlphaSites, sd = sigmaAlphaSites)
    }

    # # detection process visit effect
    # for(visit.tag in 1:n.visit){
    #   betaVisits[visit.tag] ~ dnorm(muVisits, sd = sigmaVisits)
    # }

    # Site effect
    for(site.tag in 1:n.sites){
      eta.lam[site.tag] ~ dnorm(0, sd = sigmaSites)
    }

    #Detection probability
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        #for(visit.tag in 1:n.visit){
          logit(p.tag[site.tag, spe.tag]) <- alphaSites[site.tag] + alphaSpecies[spe.tag] #+  betaVisits[visit.tag]
      #  }
      }
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
    #
    # #correlation matrix
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

  #derived quantity
  mae <- sum(abs(shan[1:n.sites] - truth[1:n.sites]))/n.species

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
                NLatent = ceiling((dim_data[2])/2),
                a1 = ceiling((dim_data[2])/2),
                a2= ceiling((dim_data[2])/2) + 5,
                a.sigma = 1,
                b.sigma = 1
  )

  pa_data <- function(data){
    zst <- apply(data, c(1,2), function(x) max(x, na.rm = TRUE))
    #zst <- data
    zst[zst!= 0] <- 1
    zst[zst == 0] <- NA
    #zst[zst!= 0] <- 1
    return(zst)
  }
  zst = pa_data(simulations_all$mat.species)

  #truth
  #if(method == "covariate" & model == "IG"){
  #  truth <- log(simulations_all$hillsGro1)
  #}else{
    truth <- log(simulations_all$hillsIDM1)
  #}

  Rmat <- diag(const$n.species) # In the paper, this is initial covariance matrix
  df <- const$n.species + 1
  #omega <- LaplacesDemon:: rwishart(df, Rmat)
  #latitude <- as.numeric(scale(simulations_all[[4]]))
  #Data for the model
  idm_data <- list(Y = simulations_all[[2]],
                   X = simulations_all[[1]],
                   detection_cov= data$detection_cov,
                   latitude = data$ecological_cov,
                   R = Rmat,
                   df = df,
                   z = zst, truth = truth)



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
                               sigmaSites = 1,
                               muSpecies = 0,
                               muLatitude = 0,
                               muVisits = 0,
                               #z= zst,
                               delta=delta,
                               phi=phi,
                               rho = 0.8,
                               mulambda = 0,
                               lamLatent= lamLatent ,
                               omega=diag(1, const$n.species),
                               sig = rgamma(const$n.species, 2,1),
                               eta.lam = rep(0, const$n.sites))

  }
  initsList <- idm_inits()


  #############################################################
  #                 Compile Model                             #
  #############################################################


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
                                         "sigmaSites",
                                         "sigmaSpecies",
                                         "sigmaLatitude",
                                         "sigmaAlphaSpecies",
                                         "sigmaAlphaSites",
                                         "sigmaLambda",
                                         "rho",
                                         "shan",
                                         "mae",
                                         #"alphaSites",
                                         "mulambda"
                                         #"lamLatent",
                                         #"phi",
                                         #"sig",
                                         #"delta",
                                         #"betaLambda",
                                         # "z"
                            ))

  # mcmcconf$removeSamplers(c("lamLatent",
  #                           "betaSpecies", "betaLatitude", "alphaSpecies", "alphaSites", "betaVisits"), print = FALSE)
  # #mcmcconf$addSampler(target = "delta", type = "RW_block", print = TRUE)
  # #mcmcconf$addSampler(target = "sig", type = "RW_block", print = TRUE)
  # #mcmcconf$addSampler(target = "phi", type = "RW_block", print = TRUE)
  # mcmcconf$addSampler(target = "lamLatent", type = "RW_block", print = TRUE)
  # mcmcconf$addSampler(target = "betaSpecies", type = "RW_block", print = TRUE)
  # mcmcconf$addSampler(target = "betaLatitude", type = "RW_block", print = TRUE)
  # mcmcconf$addSampler(target = "alphaSpecies", type = "RW_block", print = TRUE)
  # mcmcconf$addSampler(target = "alphaSites", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "betaVisits", type = "RW_block", print = TRUE)
  #mcmcconf$addSampler(target = "betaVisits", type = "AF_slice", print = TRUE)
  #build MCMC and compile
  Rmcmc <- buildMCMC(mcmcconf)

  cmcmc <- compileNimble(Rmcmc,
                         project = Cmwtc)

  #run MCMC
  start_time <- Sys.time()
  estimate <- runMCMC(cmcmc,
                      niter = 30000,#300000,
                      nburnin = 20000, #200000,
                      inits = initsList,
                      nchains = 3,
                      #thin = 20,
                      summary = TRUE,
                      samples = TRUE,
                      setSeed = TRUE,
                      samplesAsCodaMCMC = TRUE)

  end.time <- Sys.time()
  time_taken <- as.numeric(round(end.time-start_time,2))

  ################
  #Checking convergence of simulation
  # Return RHat
  ###############
  # mcmclist <- ggs(estimate$samples)
  #
  # aa <- mcmclist%>%
  #   filter(grepl('mu|sigma',Parameter))
  #
  # subset_parameters <- unique(aa$Parameter)
  #
  # # Rhat values
  # subset_Rhat <- aa%>%
  #   ggs_Rhat()
  #
  # Rhat_data <- subset_Rhat$data[,5]
  # all_rhat <- all(Rhat_data < 1.1) #Rhat is less than 1.1
  # N_over_rhat <-length(which(Rhat_data > 1.1))/length(subset_parameters) #Rhats over 1.1
  estimation <- estimate$summary$all.chains

  return(list(estimation, time_taken))
  #return(estimation)
}


# model = c("IDM", "IDM", "Spe", "IG", "IG")
# method = c("shared", "covariate", "shared", "shared", "covariate")
# ret <- list()
# for(i in 1:length(model)){
#   ret[[i]] <- run_nimble_model(simulations_all_na[[1]], method = method[i], model = model[i])
# }
