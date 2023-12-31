#This code is used to fit the MCMC to the simulated dataset



# Load packages needed to run the code
library(pbapply)
library(nimble)
library(MCMCglmm)
library(coda)
library(MCMCpack)
library(boot)
library(MASS)
library(parallel)
require(ggmcmc)

op <- pboptions(type="timer")

# Function that estimates the parameters
#Needed to run the code in parallel
run_simulations <- function(simulations_all, type, model_selection){
  library(pbapply)
  library(nimble)
  library(MCMCglmm)
  library(coda)
  library(MCMCpack)
  library(boot)
  library(MASS)
  library(parallel)
  require(ggmcmc)

  my_proportions <- function(x){
    ret <- proportions(x)
    return(ret)
  }

  nimble_proportions <- nimble::nimbleRcall(
    prototype = function(
    x=double(1)
    # beta is a vector
    ) {},
    returnType = double(1), # outcome is a vector
    Rfun = 'my_proportions'
  )

  start_time <- Sys.time()

  code <- nimbleCode({

    # Priors for observation process parameters
    for(spe.tag in 1:(n.species-1)){
      beta0[spe.tag] ~ dnorm(0,sd=10)
      beta2[spe.tag] ~ dnorm(0,sd=10)
    }

    for(spe.tag in 1:(n.species-1)){
      for (cov.tag in 1: n.cov){
        beta1[spe.tag, cov.tag] ~ dnorm(0, sd=10)
      }
    }

    #main model covariate selection
    if(model_selection == TRUE){
      p ~ dunif(0,1)
      psi ~ dbern(p)
    }else{
      p <- 1
      psi <- 1
    }


    if(type == "variable" | type == "constant" | type == "intercept"){

      for(spe.tag in 1:n.species){
        for(report.tag in 1:(n.reported-1)){
          gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
          gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }

      for(site.tag in 1:n.sites){
        lambda[n.species,site.tag] <- 1
        for(spe.tag in 1:(n.species-1)){
          #lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + inprod(beta1[spe.tag, 1:n.cov], cov[site.tag, 1:n.cov]))
          lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag, 1]*  cov[site.tag, 1] + beta1[spe.tag, 2]*  cov[site.tag, 2])
        }
      }

      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }
    }


    # Proportion for the multinomial distribution


    if(type == "variable"){

      #Prior for the Confusion Matrix
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, site.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag] + psi * gamma1[spe.tag, report.tag] * cov_omega[site.tag])
          }
        }
      }
      # Confusion Matrix for the Species
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            omega[spe.tag, report.tag, site.tag ] <- alpha[spe.tag, report.tag, site.tag ]/sum(alpha[spe.tag, 1:n.reported, site.tag ])
          }
        }
      }

      # Verified data
      for(site.tag in 1:n.sites){
        for(visit.tag in 1:n.visit){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }

      # Reported observations
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:n.reported, site.tag])
        }
      }
    }


    if(type == "intercept"){
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, site.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag])
          }
        }


      }
      # Confusion Matrix for the Species
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            omega[spe.tag, report.tag, site.tag ] <- alpha[spe.tag, report.tag, site.tag ]/sum(alpha[spe.tag, 1:n.reported, site.tag ])
          }
        }
      }

      # Verified data
      for(site.tag in 1:n.sites){
        for(visit.tag in 1:n.visit){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }

      # Reported observations
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:n.reported, site.tag])
        }
      }
    }


    if(type == "constant"){

      for(spe.tag in 1:n.species){
        for(report.tag in 1:n.reported){
          alpha[spe.tag, report.tag] ~ dexp(1)
        }
      }
      # Confusion Matrix for the Species
      for(spe.tag in 1:n.species){
        omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
      }

      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag, site.tag] ~ dcat(omega[C[visit.tag, site.tag],1:n.reported])
        }
      }
    }

    if(type == "main"){

      for(spe.tag in 1:n.species){
        for(report.tag in 1:(n.reported-1)){
          gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
          gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }


      for(site.tag in 1:n.sites){
        lambda[n.species,site.tag] <- 1
        for(spe.tag in 1:(n.species-1)){
          #lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + inprod(beta1[spe.tag, 1:n.cov], cov[site.tag, 1:n.cov]) + beta2[spe.tag]*cov_omega[site.tag])
          lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag, 1]*  cov[site.tag, 1] + beta1[spe.tag, 2]*  cov[site.tag, 2]+ beta2[spe.tag]*cov_omega[site.tag])
          }
      }

      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }

      for(spe.tag in 1:n.species){
        for(report.tag in 1:n.reported){
          alpha[spe.tag, report.tag] ~ dexp(1)
        }
      }
      # Confusion Matrix for the Species
      for(spe.tag in 1:n.species){
        omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
      }

      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag, site.tag] ~ dcat(omega[C[visit.tag, site.tag],1:n.reported])
        }
      }
    }

    if(type == "fixed-covariate"){

      for(site.tag in 1:n.sites){
        lambda[n.species,site.tag] <- 1
        for(spe.tag in 1:(n.species-1)){
          #lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + inprod(beta1[spe.tag, 1:n.cov], cov[site.tag, 1:n.cov]))
          lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag, 1]*  cov[site.tag, 1]  + beta1[spe.tag, 2]*  cov[site.tag, 2])
        }
      }

      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }


      for(spe.tag in 1:n.species){
        gamma1[spe.tag] ~ dnorm(0, sd = 1)
        for(report.tag in 1:(n.reported-1)){
          gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }


      #Prior for the Confusion Matrix
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, site.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag] + psi * gamma1[spe.tag] * cov_omega[site.tag])
          }
        }
      }

      # Confusion Matrix for the Species
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            omega[spe.tag, report.tag, site.tag ] <- alpha[spe.tag, report.tag, site.tag ]/sum(alpha[spe.tag, 1:n.reported, site.tag ])
          }
        }
      }

      # Verified data
      for(site.tag in 1:n.sites){
        for(visit.tag in 1:n.visit){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }

      # Reported observations
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:n.reported, site.tag])
        }
      }

    }


    if(type == "fixed-intercov"){

      for(site.tag in 1:n.sites){
        lambda[n.species,site.tag] <- 1
        for(spe.tag in 1:(n.species-1)){
          #lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + inprod(beta1[spe.tag, 1:n.cov], cov[site.tag, 1:n.cov]))
          lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag, 1]*  cov[site.tag, 1] + beta1[spe.tag, 2]*  cov[site.tag, 2])
        }
      }

      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }


      for(spe.tag in 1:n.species){
        gamma1[spe.tag] ~ dnorm(0, sd = 1)
        #for(report.tag in 1:(n.reported-1)){
        gamma0[spe.tag, spe.tag] <- beta2[1]

        #}
      }
      gamma0[1, 2]~ dnorm(0, sd = 1)
      gamma0[2, 1]~ dnorm(0, sd = 1)


      #Prior for the Confusion Matrix
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, site.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, site.tag ] <- exp(gamma0[spe.tag, report.tag] + psi * gamma1[spe.tag] * cov_omega[site.tag])
          }
        }
      }
      # Confusion Matrix for the Species
      for(site.tag in 1: n.sites){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            omega[spe.tag, report.tag, site.tag ] <- alpha[spe.tag, report.tag, site.tag ]/sum(alpha[spe.tag, 1:n.reported, site.tag ])
          }
        }
      }

      # Verified data
      for(site.tag in 1:n.sites){
        for(visit.tag in 1:n.visit){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }

      # Reported observations
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag,site.tag] ~ dcat(omega[C[visit.tag,site.tag],1:n.reported, site.tag])
        }
      }

    }

    if(type == "ML"){
      for(site.tag in 1:n.sites){
        lambda[n.reported,site.tag] <- 1
        for(report.tag in 1:(n.reported-1)){
          #lambda[report.tag, site.tag] <- exp(beta0[report.tag] + inprod(beta1[report.tag, 1:n.cov], cov[site.tag, 1:n.cov]))
          lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag, 1]*  cov[site.tag, 1]  + beta1[spe.tag, 2]*  cov[site.tag, 2])
        }
      }

      for(site.tag in 1:n.sites){
        for(report.tag in 1:n.reported){
          prop[report.tag,site.tag] <- lambda[report.tag, site.tag]/sum(lambda[1:n.reported, site.tag])
        }
      }

      for(site.tag in 1:n.sites){
        for(report.tag in 1:n.reported){
          prob[report.tag, site.tag] <- ((z.omega[report.tag, site.tag] + 0.01) * prop[report.tag,site.tag] )/sum((z.omega[1:n.reported, site.tag] +0.01)* prop[1:n.reported,site.tag])
        }
      }

      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          C[visit.tag,site.tag] ~ dcat(prob[1:n.reported, site.tag])
        }
      }

      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag,site.tag] ~ dcat(prop[1:n.reported, site.tag])
        }
      }
    }


    for(validation_site.tag in 1:n.validated){
      #ind_all_predictions[validation_site.tag] <- equals(Y[1, (index_for_splitting[validation_site.tag])], Validation_Y[ validation_site.tag])
      ind_all_predictions[validation_site.tag] <- equals(C[1, (index_for_splitting[validation_site.tag])], Validation_C[ validation_site.tag])
    }

    if(n.mismatch > 1){
    for(mismatch_site.tag in 1:n.mismatch){
      ind_all_precisions[mismatch_site.tag] <- equals(C[1, (validation_indices_for_mismatch[mismatch_site.tag])], validation_mismatchC[mismatch_site.tag])
    }}

    if(n.mismatch == 1){
        ind_all_precisions <- equals(C[1, (validation_indices_for_mismatch)], validation_mismatchC)
      }

    if(n.correctmatch > 1){
    for(correctmatch_site.tag in 1:n.correctmatch){
      ind_all_recall[correctmatch_site.tag] <- equals(C[1, (validation_correct_match[correctmatch_site.tag])], validation_correctmatchC[correctmatch_site.tag])
    }
    }

    if(n.correctmatch == 1){
     # for(correctmatch_site.tag in 1:n.correctmatch){
        ind_all_recall <- equals(C[1, (validation_correct_match)], validation_correctmatchC)
    #  }
    }

    accuracy <- sum(ind_all_predictions[1:n.validated])/ n.validated
    if(n.mismatch > 1){
    precision <- sum(ind_all_precisions[1:n.mismatch]) / n.mismatch
    }
    if(n.mismatch == 1){
      precision <- ind_all_precisions
    }
    if(n.mismatch == 0){
      precision <- NA
    }

    if(n.correctmatch > 1){
    recall <- sum(ind_all_recall[1:n.correctmatch])/n.correctmatch
    }

    if(n.correctmatch == 1){
  recall <- ind_all_recall
    }

    if(n.correctmatch == 0){
      recall <- NA
    }

  })

  # Retrieving data from the simulated data

  C <- simulations_all$C # Unverified citizen science observation
  Y <- simulations_all$Y #Verified Citizen science observation
  Validation_C <- simulations_all$Validation_C
  Validation_Y <- simulations_all$Validation_Y
  validation_indices_for_mismatch <- simulations_all$validation_indices_for_mismatch
  index_for_splitting <- simulations_all$index_for_splitting
  validation_mismatchC <- simulations_all$validation_mismatchC

  validation_correct_match <- simulations_all$validation_indices_for_match
  validation_correctmatchC <- simulations_all$validation_matchC
  cov_corr <- simulations_all$cov_corr
  covariate_group <- simulations_all$covariate_group

  if(covariate_group == 1){
    cov <- t(simulations_all$cov) #Covariates
    cov_omega <- simulations_all$cov_omega
  }else{
    cov <- t(rbind(simulations_all$cov[1,], cov_corr[,1])) #Covariates
    cov_omega <- cov_corr[,2]
  }


  if(type == "ML"){
    p.omega <- t(simulations_all$omega)
  }else{
    p.omega <- NULL
  }

  data <- list(C,Y,
               cov,cov_omega,
               Validation_C,Validation_Y,
               validation_indices_for_mismatch, index_for_splitting,
               validation_mismatchC,  validation_correct_match,
               validation_correctmatchC, p.omega)

  #Constants for the model
  dim_data <- dim(data[[1]])
  n.cov = ncol(cov)
  const <- list(n.sites=dim_data[2],
                n.species = 2, #max(c(data[[1]]), na.rm = TRUE), #fix them for consistency
                n.reported = 3, #max(c(data[[2]]), na.rm = TRUE),
                n.visit=(dim_data[1]),
                n.cov=n.cov,
                #nref = max(c(data[[1]]))-1,
                n.mismatch = length(data[[7]]),
                n.validated = length(data[[8]]),
                n.correctmatch = length(data[[10]]),
                validation_indices_for_mismatch = data[[7]],
                index_for_splitting = data[[8]],
                validation_correct_match = data[[10]],
                z.omega = data[[12]]
  )

  # Data for the Model

  idm_data <- list(C = data[[1]],
                   Y = data[[2]],
                   cov = data[[3]],
                   cov_omega = as.numeric(data[[4]]),
                   Validation_C= data[[5]],
                   # Validation_Y = data[[6]],
                   validation_mismatchC = data[[9]],
                   validation_correctmatchC = data[[11]])

  # Initials for the model
  if(type == "variable" | type == "intercept"){
    alpha <- omega <- array(NA, dim = c(const$n.species, const$n.reported, const$n.sites))

    for(site.tag in 1:const$n.sites){
      for(spe.tag in 1:const$n.species){
        for(report.tag in 1:(const$n.reported)){
          alpha[spe.tag,report.tag, site.tag] <- 1
        }
      }
    }

    for(site.tag in 1:(const$n.sites)){
      for(spe.tag in 1:(const$n.species)){
        omega[spe.tag,, site.tag]<- proportions(alpha[spe.tag,1:(const$n.reported),site.tag])
      }
    }

    # idm_inits <- function(){list(beta1 = matrix(rnorm(const$n.cov * (const$n.species -1)), ncol = const$n.cov),
    #                              omega=omega,
    #                              beta0 = rnorm((const$n.species-1)),
    #                              gamma0 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported - 1)),
    #                              gamma1 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported -1)),
    #                              alpha = alpha
    # )
    # }

    idm_inits <- function(){list(beta1 = matrix(1, nrow= (const$n.species -1), ncol = const$n.cov),
                                 omega=omega,
                                 beta0 = rep(1, (const$n.species-1)),
                                 gamma0 = matrix(1, nrow= const$n.species, ncol=(const$n.reported - 1)),
                                 gamma1 = matrix(1, nrow= const$n.species, ncol=(const$n.reported -1)),
                                 alpha = alpha
    )
    }
  }else{if(type =="fixed-intercov"| type == "fixed-covariate"){
    alpha <- omega <- array(NA, dim = c(const$n.species, const$n.reported, const$n.sites))

    for(site.tag in 1:const$n.sites){
      for(spe.tag in 1:const$n.species){
        for(report.tag in 1:(const$n.reported)){
          alpha[spe.tag,report.tag, site.tag] <- 1
        }
      }
    }

    for(site.tag in 1:(const$n.sites)){
      for(spe.tag in 1:(const$n.species)){
        omega[spe.tag,, site.tag]<- proportions(alpha[spe.tag,1:(const$n.reported),site.tag])
      }
    }

    idm_inits <- function(){list(beta1 = matrix(rnorm(const$n.cov * (const$n.species -1)), ncol = const$n.cov),
                                 omega=omega,
                                 beta0 = rnorm((const$n.species-1)),
                                 gamma0 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported - 1)),
                                 gamma1 = rnorm(const$n.species),
                                 alpha = alpha
    )
    }
  }else{if(type == "ML"){
    idm_inits <- function(){list(beta1 = matrix(1, nrow = (const$n.reported -1), ncol = const$n.cov),
                                 #z.omega=z.omega,
                                 beta0 = rep(1, (const$n.reported-1))
    )
    }}else{
      alpha <- matrix(NA, nrow=const$n.species, ncol=(const$n.reported))
      for(spe.tag in 1:const$n.species){
        for(k in 1:(const$n.reported)){
          alpha[spe.tag,k] <- 1
        }
      }


      #alpha <- rep(1, (const$n.species+1))


      # Initials for the model
      omega <- matrix(NA, nrow=const$n.species,
                      ncol = (const$n.reported))
      for(spe.tag in 1:(const$n.species)){
        omega[spe.tag,]<- extraDistr::rdirichlet(1, alpha = alpha[spe.tag,1:(const$n.reported)])
      }

      idm_inits <- function(){list(beta1 = matrix(rnorm(const$n.cov * (const$n.species -1)), ncol = const$n.cov),
                                   omega=omega,
                                   beta0 = rep(1, (const$n.species-1)),
                                   beta2 = rep(1, (const$n.species-1)),
                                   gamma0 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported - 1)),
                                   gamma1 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported -1)),
                                   alpha = alpha
      )
      }
    }
  }
  }
  #}
  initsList <- idm_inits()

  #Putting all together for the creating compilation
  # modelInfo <- list(
  #   code = code,
  #   constants = const,
  #   data = idm_data,
  #   inits = initsList
  # )

  #Create the model in nimble
  message("Setting up the NIMBLE model")
  mwtc <- nimble::nimbleModel(code,
                              data = idm_data,
                              constants = const,
                              inits = initsList)

  # Create the model in C
  message("Create the model in C")
  Cmwtc <- nimble::compileNimble(mwtc,
                                 showCompilerOutput = FALSE)


  message("MCMC configuration")
  mcmcconf <- nimble::configureMCMC(Cmwtc,
                                    monitors = c("beta0","beta1",
                                                 "gamma0","gamma1",
                                                 "accuracy", "precision",
                                                 "p", "recall"),
                                    enableWAIC = TRUE)
  message("Builing the MCMC model")
  Rmcmc <- nimble::buildMCMC(mcmcconf)

  # Compile
  message("Compiling the model")
  cmcmc <- nimble::compileNimble(Rmcmc,
                                 project = Cmwtc)

  # Run the MCMC
  message("Running the MCMC")
  #ret <- cmcmc$run(niter = 1000)
  # Run the MCMC
  mcmc.out <- runMCMC(cmcmc,
                      niter = 10000,
                      nchains = 3,
                      nburnin = 5000,
                      inits = initsList,
                      #thin = 5,
                      setSeed = TRUE,
                      samples=TRUE,
                      samplesAsCodaMCMC = TRUE,
                      summary = TRUE,
                      WAIC = TRUE)

  # Return the summary of MCMC results
 output <-  mcmc.out$summary$all.chains
  end.time <- Sys.time()
  time_taken <- as.numeric(round(end.time-start_time,2))
  iteration <- simulations_all$iteration
  waic = mcmc.out$WAIC

  return(list(output = output,
              n.mismatch = const$n.mismatch,
              n.validated = const$n.validated,
              n.correctmatch = const$n.correctmatch,
              waic = waic,
              iteration = iteration))
}


