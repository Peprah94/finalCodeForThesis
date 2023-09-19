#This script fits the model using MCMC. 

#load the data
load("nimble_data_ML.RData")

# Required packages
library(spatstat)
library(dplyr)
library(sp)
library(gstat)
library(ggplot2)
library(INLA)
library(rgeos)
library(pbapply)
library(nimble)
library(MCMCglmm)
library(coda)
library(MCMCpack)
library(boot)
library(MASS)
library(parallel)
require(ggmcmc)

run_simulations <- function(simulations_all, type, model_selection){
  # packages copied within the code for the parallelisation
  library(nimble)
  library(MCMCglmm)
  library(coda)
  library(mcmcplots)
  library(MCMCpack)
  library(boot)
  #library(extraDistr)
  library(MASS)
  #library(spdep)
  require(ggmcmc)
  
  nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
  
  my_proportions <- function(x){
    ret <- proportions(x)
    return(ret)
  }
  
  nimble_proportions <- nimbleRcall(
    prototype = function(
    x=double(1)
    # beta is a vector
    ) {},
    returnType = double(1), # outcome is a vector
    Rfun = 'my_proportions'
  )
  
  start_time <- Sys.time()
  
  
  
  code <- nimbleCode({
    # Priors for true state observation process
    for(spe.tag in 1:(n.species-1)){
      beta0[spe.tag] ~ dnorm(0,sd=10) #intercept
      beta5[spe.tag] ~ dnorm(0, sd=10) #effect of covariate in main model
    }
    
    for(spe.tag in 1:(n.species-1)){
      for (cov.tag in 1: n.cov){
        beta1[spe.tag, cov.tag] ~ dnorm(0, sd=10) #covariate effect
      }
    }
    
    # Bayesian variable selection probability
    if(model_selection == TRUE){
      p ~ dunif(0,1)
      psi ~ dbern(p)
    }else{
      p <- 1
      psi <- 1
    }
    
    #Classification process prior
    if(type == "variable" | type == "constant" | type == "intercept" ){
      
      # classification process covariates
      #Note gamma is omega in the main paper. 
      #We use gamma to avoid confusions in coding
      for(spe.tag in 1:n.species){
        for(report.tag in 1:(n.reported-1)){
          gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
          gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }
      
      #True state observation process intensity
      if(n.cov == 1){
      for(site.tag in 1:n.sites){
        lambda[n.species,site.tag] <- 1
        for(spe.tag in 1:(n.species-1)){
          lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag, 1]* cov[1, site.tag])
        }
      }
      }else{
        for(site.tag in 1:n.sites){
          lambda[n.species,site.tag] <- 1
          for(spe.tag in 1:(n.species-1)){
            lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + inprod(beta1[spe.tag, 1:n.cov], cov[1:n.cov, site.tag]))
          }
        }
      }
      
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }
    }
    
  
    ##############################
    # Specifying classification process model for all the study scenarios in Table 2
    ################################
    if(type == "variable"){
      
      #Prior for the Confusion Matrix
      for(user.tag in 1: n.users){
        user_effect[user.tag] ~ dnorm(0, tau = 0.001)
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, user.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
           alpha[spe.tag, report.tag, user.tag ] <- exp(gamma0[spe.tag, report.tag] + psi * gamma1[spe.tag, report.tag] * users_id[user.tag])
            #alpha[spe.tag, report.tag, user.tag ] <- exp(gamma0[spe.tag, report.tag] + psi *  user_effect[users_id[user.tag]])
          }
        }
      }
      # Confusion Matrix for the Species
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            z.omega[spe.tag, report.tag, user.tag ] <- alpha[spe.tag, report.tag, user.tag ]/sum(alpha[spe.tag, 1:n.reported, user.tag ])
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
          Y[visit.tag,site.tag] ~ dcat(z.omega[C[visit.tag,site.tag],1:n.reported, users[site.tag] ])
        }
      }
    }
    
    
    if(type == "intercept"){
      #Prior for the Confusion Matrix
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, user.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
            alpha[spe.tag, report.tag, user.tag ] <- exp(gamma0[spe.tag, report.tag] )
          }
        }
      }
      # Confusion Matrix for the Species
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            z.omega[spe.tag, report.tag, user.tag ] <- alpha[spe.tag, report.tag, user.tag ]/sum(alpha[spe.tag, 1:n.reported, user.tag ])
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
          Y[visit.tag,site.tag] ~ dcat(z.omega[C[visit.tag,site.tag],1:n.reported, users[site.tag] ])
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
        z.omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
      }
      
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag, site.tag] ~ dcat(z.omega[C[visit.tag, site.tag],1:n.reported])
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
      
      for(user.tag in 1:n.users){
        user_effect[user.tag] ~ dnorm(0, sd = 1)
      }
      
      if(n.cov == 1){
        for(site.tag in 1:n.sites){
          lambda[n.species,site.tag] <- 1
          for(spe.tag in 1:(n.species-1)){
            lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag, 1]* cov[1, site.tag] + psi*user_effect[users[site.tag]])
          }
        }
      }else{
        for(site.tag in 1:n.sites){
          lambda[n.species,site.tag] <- 1
          for(spe.tag in 1:(n.species-1)){
            lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + inprod(beta1[spe.tag, 1:n.cov], cov[1:n.cov, site.tag]) + psi*user_effect[users[site.tag]])
          }
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
        z.omega[spe.tag, 1:n.reported] ~ ddirch(alpha[spe.tag,1:n.reported])
      }
      
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          C[visit.tag, site.tag] ~ dcat(prop[1:n.species,site.tag])
        }
      }
      for(visit.tag in 1:n.visit){
        for(site.tag in 1:n.sites){
          Y[visit.tag, site.tag] ~ dcat(z.omega[C[visit.tag, site.tag],1:n.reported])
        }
      }
    }
    
    if(type == "fixed-covariate"){
      
      for(spe.tag in 1:n.species){
        for(report.tag in 1:(n.reported-1)){
          gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
          gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }
      
      if(n.cov == 1){
        for(site.tag in 1:n.sites){
          lambda[n.species,site.tag] <- 1
          for(spe.tag in 1:(n.species-1)){
            lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag, 1]* cov[1, site.tag])
          }
        }
      }else{
        for(site.tag in 1:n.sites){
          lambda[n.species,site.tag] <- 1
          for(spe.tag in 1:(n.species-1)){
            lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + inprod(beta1[spe.tag, 1:n.cov], cov[1:n.cov, site.tag]))
          }
        }
      }
      
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }
      
      for(user.tag in 1: n.users){
        user_effect[user.tag] ~ dnorm(0, tau = 0.001)
        for(spe.tag in 1:n.species){
          # for(user.tag in 1:n.users){
          alpha[spe.tag, n.reported, user.tag] <- 1
          for(report.tag in 1:(n.reported-1)){
           # alpha[spe.tag, report.tag, user.tag] <- exp(gamma0[spe.tag, report.tag] + psi *  user_effect[users_id[user.tag]])
            alpha[spe.tag, report.tag, user.tag ] <- exp(gamma0[spe.tag, report.tag] + psi * gamma1[spe.tag, spe.tag] * users_id[user.tag])
            # }
          }
        }
      }
      # Confusion Matrix for the Species
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            #for(user.tag in 1:n.users){
            z.omega[spe.tag, report.tag, user.tag] <- alpha[spe.tag, report.tag, user.tag]/sum(alpha[spe.tag, 1:n.reported, user.tag])
            #  }
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
          Y[visit.tag,site.tag] ~ dcat(z.omega[C[visit.tag,site.tag],1:n.reported, users[site.tag] ])
        }
      }
    }
    
    
    
    if(type == "fixed-intercov"){
      
      if(n.cov == 1){
        for(site.tag in 1:n.sites){
          lambda[n.species,site.tag] <- 1
          for(spe.tag in 1:(n.species-1)){
            lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + beta1[spe.tag, 1]* cov[1, site.tag])
          }
        }
      }else{
        for(site.tag in 1:n.sites){
          lambda[n.species,site.tag] <- 1
          for(spe.tag in 1:(n.species-1)){
            lambda[spe.tag, site.tag] <- exp(beta0[spe.tag] + inprod(beta1[spe.tag, 1:n.cov], cov[1:n.cov, site.tag]))
          }
        }
      }
      
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          prop[spe.tag,site.tag] <- lambda[spe.tag, site.tag]/sum(lambda[1:n.species, site.tag])
        }
      }
      
      for(spe.tag in 1:n.species){
        gamma0[spe.tag, spe.tag] <- beta5[1]
      }
      gamma0[1, 2]~ dnorm(0, sd = 1)
      gamma0[1, 3]~ dnorm(0, sd = 1)
      #gamma0[1, 4]~ dnorm(0, sd = 1)
     # gamma0[1, 5]~ dnorm(0, sd = 1)
      
      gamma0[2, 1]~ dnorm(0, sd = 1)
      gamma0[2, 3]~ dnorm(0, sd = 1)
      #gamma0[2, 4]~ dnorm(0, sd = 1)
      #gamma0[2, 5]~ dnorm(0, sd = 1)
      
      gamma0[3, 1]~ dnorm(0, sd = 1)
      gamma0[3, 2]~ dnorm(0, sd = 1)
      #gamma0[3, 4]~ dnorm(0, sd = 1)
     # gamma0[3, 5]~ dnorm(0, sd = 1)
      
      #gamma0[1, 2]~ dnorm(0, sd = 10)
      #gamma0[1, 2]~ dnorm(0, sd = 10)
      
      for(spe.tag in 1:n.species){
        for(report.tag in 1:(n.reported-1)){
          #gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
          gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }
      
      #Prior for the Confusion Matrix
      for(user.tag in 1: n.users){
        user_effect[user.tag] ~ dnorm(0, tau = 0.001)
        for(spe.tag in 1:n.species){
          alpha[spe.tag, n.reported, user.tag ] <- 1
          for(report.tag in 1:(n.reported-1)){
           # alpha[spe.tag, report.tag, user.tag] <- exp(gamma0[spe.tag, report.tag] + psi *  user_effect[users_id[user.tag]])
            alpha[spe.tag, report.tag, user.tag ] <- exp(gamma0[spe.tag, report.tag] + psi * gamma1[spe.tag, report.tag] * users_id[user.tag])
          }
        }
      }
      
      # Confusion Matrix for the Species
      for(user.tag in 1: n.users){
        for(spe.tag in 1:n.species){
          for(report.tag in 1:n.reported){
            z.omega[spe.tag, report.tag, user.tag ] <- alpha[spe.tag, report.tag, user.tag ]/sum(alpha[spe.tag, 1:n.reported, user.tag ])
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
          Y[visit.tag,site.tag] ~ dcat(z.omega[C[visit.tag,site.tag],1:n.reported, users[site.tag]])
        }
      }
      
    } 
    
    if(type == "ML"){
      for(report.tag in 1:(n.reported-1)){
        for (cov.tag in 1: n.cov){
          beta2[report.tag, cov.tag] ~ dnorm(0, sd=10)
        }
      }
      
      if(n.cov == 1){
        for(site.tag in 1:n.sites){
          lambda[n.reported,site.tag] <- 1
          for(report.tag in 1:(n.reported-1)){
            lambda[report.tag, site.tag] <- exp(beta0[report.tag] + beta2[report.tag, 1]* cov[1, site.tag])
          }
        }
      }else{
        for(site.tag in 1:n.sites){
          lambda[n.reported,site.tag] <- 1
          for(report.tag in 1:(n.reported-1)){
            lambda[report.tag, site.tag] <- exp(beta0[report.tag] + inprod(beta2[report.tag, 1:n.cov], cov[1:n.cov, site.tag]))
          }
        }
      }

      
      z.omega <- 1
      for(site.tag in 1:n.sites){
        for(report.tag in 1:n.reported){
          prop[report.tag,site.tag] <- lambda[report.tag, site.tag]/sum(lambda[1:n.reported, site.tag])
        }
      }
      
      for(spe.tag in 1:n.species){
        for(report.tag in 1:(n.reported-1)){
          gamma0[spe.tag, report.tag] ~ dnorm(0, sd = 1)
          gamma1[spe.tag, report.tag] ~ dnorm(0, sd = 1)
        }
      }
      
      for(site.tag in 1:n.sites){
        for(report.tag in 1:n.reported){
          prob[report.tag, site.tag] <- ((fixed.omega[report.tag, site.tag] + 0.01) * prop[report.tag,site.tag] )/sum((fixed.omega[1:n.reported, site.tag] +0.01)* prop[1:n.reported,site.tag])
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
      ind_all_predictions[validation_site.tag] <- equals(C[1, (index_for_splitting[validation_site.tag])], Validation_C[ validation_site.tag])
    }
    
    for(mismatch_site.tag in 1:n.mismatch){
      ind_all_precisions[mismatch_site.tag] <- equals(C[1, (validation_indices_for_mismatch[mismatch_site.tag])], validation_mismatchC[mismatch_site.tag])
    }
    
    for(correctmatch_site.tag in 1:n.correctmatch){
      ind_all_recall[correctmatch_site.tag] <- equals(C[1, (validation_correct_match[correctmatch_site.tag])], validation_correctmatchC[correctmatch_site.tag])
    }
    
    # returning prediction metrics
    accuracy <- sum(ind_all_predictions[1:n.validated])/ n.validated
    precision <- sum(ind_all_precisions[1:n.mismatch]) / n.mismatch
    recall <- sum(ind_all_recall[1:n.correctmatch])/n.correctmatch
  })
  
  C <- simulations_all$C # Unverified citizen science observation
  Y <- simulations_all$Y #Verified Citizen science observation
  Validation_C <- simulations_all$Validation_C
  Validation_Y <- simulations_all$Validation_Y
  validation_indices_for_mismatch <- simulations_all$validation_indices_for_mismatch
  index_for_splitting <- simulations_all$index_for_splitting
  validation_mismatchC <- simulations_all$validation_mismatchC
  cov <- t(as.matrix(simulations_all$cov[2,])) #Covariates
  cov_omega <- as.numeric(simulations_all$cov_omega)
  validation_correct_match <- simulations_all$validation_indices_for_match
  validation_correctmatchC <- simulations_all$validation_matchC
  cov_corr <- simulations_all$cov_corr
  covariate_group <- simulations_all$covariate_group
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

  
  dim_data <- dim(data[[1]])
  n.cov = nrow(cov)
  const <- list(n.sites=dim_data[2],
                n.species=max(c(data[[1]]), na.rm = TRUE),
                n.reported = max(c(data[[2]]), na.rm = TRUE),
                n.visit=(dim_data[1]),
                n.cov=n.cov,
                n.users = length(unique(simulations_all$cov_omega)),
                n.mismatch = length(data[[7]]),
                n.validated = length(data[[8]]),
                n.correctmatch = length(data[[10]]),
                validation_indices_for_mismatch = data[[7]],
                index_for_splitting = data[[8]],
                validation_correct_match = data[[10]],
                fixed.omega = data[[12]],
                users_id = unique(data[[4]]),
                users = data[[4]]
  )
  
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
    alpha <- omega <- array(NA, dim = c(const$n.species, const$n.reported, const$n.user))
    
    for(user.tag in 1:const$n.users){
      for(spe.tag in 1:const$n.species){
        for(report.tag in 1:(const$n.reported)){
          alpha[spe.tag,report.tag, user.tag] <- 1
        }
      }
    }
    
    for(user.tag in 1:(const$n.users)){
      for(spe.tag in 1:(const$n.species)){
        omega[spe.tag,, user.tag]<- proportions(alpha[spe.tag,1:(const$n.reported),user.tag])
      }
    }
    

    
    idm_inits <- function(){list(beta1 = matrix(1, nrow= (const$n.species -1), ncol = const$n.cov),
                                 z.omega=omega,
                                 beta0 = rep(1, (const$n.species-1)),
                                 gamma0 = matrix(1, nrow= const$n.species, ncol=(const$n.reported - 1)),
                                 gamma1 = matrix(1, nrow= const$n.species, ncol=(const$n.reported -1)),
                                 alpha = alpha
    )
    }
  }else{if(type =="fixed-intercov"| type == "fixed-covariate"){
    alpha <- omega <- array(NA, dim = c(const$n.species, const$n.reported, const$n.users))
    
    for(user.tag in 1:const$n.users){
      for(spe.tag in 1:const$n.species){
        for(report.tag in 1:(const$n.reported)){
          alpha[spe.tag,report.tag, user.tag] <- 1
        }
      }
    }
    
    for(user.tag in 1:(const$n.users)){
      for(spe.tag in 1:(const$n.species)){
        omega[spe.tag,, user.tag]<- proportions(alpha[spe.tag,1:(const$n.reported),user.tag])
      }
    }
    
    idm_inits <- function(){list(beta1 = matrix(rnorm(const$n.cov * (const$n.species -1)), ncol = const$n.cov),
                                 z.omega=omega,
                                 beta0 = rnorm((const$n.species-1)),
                                 gamma0 = matrix(rnorm(const$n.species*(const$n.reported-1)), nrow= const$n.species, ncol=(const$n.reported - 1)),
                                 gamma1 = matrix(1, nrow= const$n.species, ncol=(const$n.reported -1)),
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
  modelInfo <- list(
    code = code,
    constants = const,
    data = idm_data,
    inits = initsList
  )
  
  #Create the model in nimble
  mwtc <- nimbleModel(code,data = idm_data, 
                      constants = const, 
                      inits = initsList)
  
  # Create the model in C
  Cmwtc <- compileNimble(mwtc, 
                         showCompilerOutput = FALSE)
  
  mcmcconf <- configureMCMC(Cmwtc, 
                            monitors = c("beta0",
                                         "beta1", "beta5",
                                         "gamma0", "gamma1",
                                         "accuracy", "precision", "p", 
                                          "recall", "z.omega"),
                            enableWAIC = TRUE)

  Rmcmc <- buildMCMC(mcmcconf)
  
  # Compile 
  cmcmc <- compileNimble(Rmcmc, 
                         project = Cmwtc,
                         resetFunctions = TRUE)

  # Run the MCMC
  mcmc.out <- runMCMC(cmcmc, 
                      niter = 10000,
                      nchains =3,
                      nburnin = 5000,
                      inits = initsList,
                      thin = 5, 
                      setSeed = TRUE, 
                      samples=TRUE, 
                      samplesAsCodaMCMC = TRUE, 
                      summary = TRUE, 
                      WAIC = TRUE)
  
  # Return the summary of MCMC results
  output <- mcmc.out$summary$all.chains
  waic <- mcmc.out$WAIC
  end.time <- Sys.time()
  time_taken <- as.numeric(round(end.time-start_time,2))
  
  mcmclist <- ggs(mcmc.out$samples)
  aa <- mcmclist%>%
    filter(grepl('gamma|beta',Parameter))
  
  subset_parameters<- unique(aa$Parameter)
  
  # Rhat values
  subset_Rhat <- aa%>%
    ggs_Rhat()
  Rhat_data <- subset_Rhat$data[,5]
  all_rhat <- all(Rhat_data < 1.1) #Rhat is less than 1.05
  N_over_rhat <-length(which(Rhat_data > 1.1))/length(subset_parameters) #Rhats over 1.04

  return(list(output = output,  N_over_rhat =  N_over_rhat, waic = waic ))
}