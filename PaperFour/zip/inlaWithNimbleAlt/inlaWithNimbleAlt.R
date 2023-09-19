#Zero Inflated Poisson
library(INLA)
library(pscl)
library(mvtnorm)
library(readr)
library(dplyr)
library(myphdthesis)

# Load data
load("zip/zinb.RData")

summary(zinb)

#Maximum likehood estimates
# Full data
d <- zinb

# ZIP
ml.res <- summary(zeroinfl(count ~ child + camper | persons, data = d))
ml.res

#estimated$cous of gamma to use as proposal
muGammaEst <- as.numeric(ml.res$coefficients[[2]][1:2,1])
covGammaEst <- diag(3*as.numeric(ml.res$coefficients[[2]][1:2,2]))

#values for use
x <- cbind(zinb$child,zinb$camper, zinb$persons)

# zinb%>%
# dplyr::select(child, camper, persons)

y <- zinb$count

# Define conditional model to be fitted using INLA
fit.inla <- function(x ,
                     y ,
                     beta) {

  data <- data.frame(count = y,
                     child = x[,1],
                     camper = x[,2],
                     persons = x[, 3])

  logit_pi <- beta[1] + beta[2] * data$persons

  # Define pi hyper for likelihood
  hyper_pi <- lapply(logit_pi, function(X) {
    list(hyper = list(prob = list(fixed = TRUE, initial = X)))
  })

  # Define COUNT as diagonal matrix of observed counts
  COUNT <- matrix(NA, nrow = nrow(data), ncol = nrow(data))
  diag(COUNT) <- data$count

  res <- inla(COUNT ~ child + camper,
              data = list(COUNT = COUNT,
                          child = data$child,
                          camper = data$camper),
              family = rep("zeroinflatedpoisson1", nrow(data)),
              num.threads = "1:1",
              control.fixed = list(prec.intercept = 0.001),
              control.family = hyper_pi,
              control.compute = list(config = TRUE),
              control.predictor = list(compute = TRUE, link = 1)
  )

  samples <- inla.posterior.sample(1, res)
  fitted_values <- samples[[1]]$latent[grepl("Predictor",rownames(samples[[1]]$latent) ),1]
  intercept = samples[[1]]$latent[grepl("Intercept",rownames(samples[[1]]$latent) ),1]
  beta1 = samples[[1]]$latent[grepl("child",rownames(samples[[1]]$latent) ),1]
  beta2 = samples[[1]]$latent[grepl("camper",rownames(samples[[1]]$latent) ),1]

  ret <- cbind(fitted_values, intercept, beta1, beta2)
  #colnames(ret) <- c("mld", fixedVals)
  return(ret)
}


nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1)
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)

 # nimbleINLA(x = x,
 #            y = y,
 #            beta = beta)

### this whole section is available at
## https://r-nimble.org/nimbleExamples/zero_inflated_poisson.html
## with documentation
##


dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), zeroProb = double(),
                 log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if(x != 0) {
      ## return the log probability if log = TRUE
      if(log) return(dpois(x, lambda, log = TRUE) + log(1-zeroProb))
      ## or the probability if log = FALSE
      else return((1-zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1-zeroProb) * dpois(0, lambda, log = FALSE)
    if(log) return(log(totalProbZero))
    return(totalProbZero)
  })

rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if(isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

registerDistributions(list(
  dZIP = list(
    BUGSdist = "dZIP(lambda, zeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'lambda = double()', 'zeroProb = double()')
  )))

code <- nimbleCode({

  #for( i in 1:3){
  #  beta[i] ~ dnorm(0, 0.001)
  #}

  gamma[1:2] ~ dmnorm(muGammaEst[1:2], cov = covGammaEst[1:2, 1:2])
  #alpha ~ dnorm(0, 0.001)
  for(i in 1:N){
    logit(p[i]) <- gamma[1] + gamma[2]* x[i,3]
  }


log(mu[1:N]) <- linPred[1:N, 1]
beta0 <- linPred[1, 2]
beta1 <- linPred[1, 3]
beta2 <- linPred[1, 4]


linPred[1:N, 1:4] <- nimbleINLA(x[1:N, 1:3], yObs[1:N], gamma[1:2])
  #Fitting the inla with the simulated parameters


  # linear model specification
  for(i in 1:N){
    y[i] ~ dZIP(mu[i], zeroProb = p[i])
  }

})

inla_data <- list(y = as.numeric(zinb$count),
                  yObs = as.numeric(zinb$count),
                  x = x)

#Constants
const <- list(N = nrow(inla_data$x),
              muGammaEst = muGammaEst,
              covGammaEst = covGammaEst)

# Initial values
idm_inits <- function(){list(
  gamma = c(1,-1)
)
}

# Fit the INLA within Nimble model
samplers <- c('RW_block')
zipModelAlt <- list()
for(i in seq_along(samplers)){
  zipModelAlt[[i]] = INLAWiNim(data = inla_data,
                                    code = code,
                                    modelData = inla_data,
                                    modelConstants = const,
                                    modelInits = idm_inits,
                                    fam = "zeroinflatedpoisson1",
                                    mcmcSamplerChange = TRUE,
                                    parametersForSamplerChange = "gamma",
                                    parametersToMonitor = c("gamma", "beta0", "beta1", "beta2"),
                                    newSampler = samplers[i],#'AF_slice',
                                    newSamplerControl = list(propCov = covGammaEst,
                                                             adaptive = FALSE),
                                    mcmcConfiguration =  list(n.chains = 1,
                                                              n.iterations = 100500,
                                                              n.burnin = 50500,
                                                              n.thin = 5,
                                                              setSeed = TRUE,
                                                              samples=TRUE,
                                                              samplesAsCodaMCMC = TRUE,
                                                              summary = TRUE,
                                                              WAIC = FALSE))
}

save(zipModelAlt, file = "zip/inlaWithNimbleAlt/zipModelAlt.RData")
