#Zero Inflated Poisson

library(INLA)
library(pscl)
library(mvtnorm)
library(readr)
library(dplyr)
library(myphdthesis)

# Load the required data
load("zip/zinb.RData")
summary(zinb)

#Maximum likehood estimates
# Full data
d <- zinb

# ZIP
ml.res <- summary(zeroinfl(count ~ child + camper | persons, data = d))
ml.res

#estimates of gamma to use as proposal
muGammaEst <- as.numeric(ml.res$coefficients[[2]][1:2,1])
covGammaEst <- diag(3*as.numeric(ml.res$coefficients[[2]][1:2,2]))

#values for use
x <- cbind(zinb$child,zinb$camper, zinb$persons)

y <- zinb$count

 family = "zeroinflatedpoisson1"
 fixedVals <- c("intercept", "beta1", "beta2")
 
 # Define conditional model to be fitted using R-INLA
fit.inla <- function(x ,
                     y ,
                     beta,
                     fixedVals,
                     #interInModel,
                     family) {
  data <- data.frame(count = y,
                     child = x[,1],
                     camper = x[,2],
                     persons = x[, 3])

  logit_pi <- beta[1] + beta[2] * data$persons

  # Define pi hyper for likelihood
  hyper_pi <- lapply(logit_pi, function(X) {
    list(hyper = list(prob = list(fixed = TRUE,
                                  initial = X)))
  })

  # Define COUNT as diagonal matrix of observed counts
  COUNT <- matrix(NA, nrow = nrow(data), ncol = nrow(data))
  diag(COUNT) <- data$count

  errorIndicator <- inherits(try(res <- inla(COUNT ~ child + camper,
              data = list(COUNT = COUNT,
                          child = data$child,
                          camper = data$camper),
              family = rep(family, nrow(data)),
              num.threads = "1:1",
              control.fixed = list(prec.intercept = 0.001),
              control.family = hyper_pi
  ), silent = TRUE),
  "try-error")
  if(errorIndicator){
    intercept = NA
    beta1 = NA
    beta2 = NA
    fitted_values = -Inf
  }else{
    
    samples <- inla.posterior.sample(1, res)
    intercept = samples[[1]]$latent[grepl("Intercept",rownames(samples[[1]]$latent) ),1]
    beta1 = samples[[1]]$latent[grepl("child",rownames(samples[[1]]$latent) ),1]
    beta2 = samples[[1]]$latent[grepl("camper",rownames(samples[[1]]$latent) ),1]
  # intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  # beta1 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[2]])
  # beta2 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[3]])
  fitted_values = res$mlik[[1]]
  }
  #fitted_values = res$summary.fitted.values[,"mean"]
  ret <- cbind(fitted_values, intercept, beta1, beta2)
  colnames(ret) <- c("mld", fixedVals)
  return(ret)
}


nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1), # beta is a vector
    fixedVals = character(1, default = c("intercept", "beta1", "beta2")),
    #interInModel = double(0, default = 1),
    family = character(0, default = "zeroinflatedpoisson1")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)

#test if the function works
# nimbleINLA(x = x,
#            y = y,
#            beta = beta)


code <- nimbleCode({

  for( i in 1:3){
    beta[i] ~ dnorm(0, 0.001)
  }

  #for( i in 1:2){
    gamma[1:2] ~ dmnorm(muGammaEst[1:2], cov = covGammaEst[1:2, 1:2])
  #}
  #alpha ~ dnorm(0, 0.001)
  for(i in 1:N){
 logit(p[i]) <- gamma[1] + gamma[2]* x[i,3]
  }


  for(i in 1:N){
    log(mu[i]) <- beta[1] + beta[2]*x[i,1] + beta[3]*x[i,2]
  }
  #Fitting the inla with the simulated parameters


  # linear model specification
  for(i in 1:N){
    y[i] ~ dZIP(mu[i], zeroProb = p[i])
  }

})

inla_data <- list(y = as.numeric(y),
                  x = as.matrix(as.data.frame(x))
)

#Constants
const <- list(N = nrow(inla_data$x),
              muGammaEst = muGammaEst,
              covGammaEst = covGammaEst)

# Initial values
idm_inits <- function(){list(beta = rep(0, 3),
                             gamma = c(1,-1)
)
}

# Define samplers
samplers <- c("RW_INLA_block", "AFSS_INLA_block", "RW_block")
inlaMCMCtype <- c("inlamcmc", "inlamcmc", "mcmc")

zipModel <- list()

for(i in seq_along(samplers)){
  zipModel <- INLAWiNimDataGenerating(data = c("y"),
                               covariate = x,
                               code = code,
                               family = "zeroinflatedpoisson1",
                               modelData = inla_data,
                               modelConstants = const,
                               modelInits = idm_inits,
                               nimbleINLA = nimbleINLA,
                               inlaMCMC = inlaMCMCtype[i],
                               inlaMCsampler = samplers[i],
                               samplerControl = list(adaptive = FALSE),
                               parametersToMonitor = list(mcmc = c("gamma"),
                                                          inla = c("beta")),
                               mcmcConfiguration = list(n.chains = 1,
                                                        n.iterations = 100500,
                                                        n.burnin = 50500,
                                                        n.thin = 5,
                                                        setSeed = TRUE,
                                                        samples=TRUE,
                                                        samplesAsCodaMCMC = TRUE,
                                                        summary = TRUE,
                                                        WAIC = FALSE)
)
  save(zipModel, file = paste0("zip/inlaWithNimble/inlaNimZIP",i,".RData"))
}




