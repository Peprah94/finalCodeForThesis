# Fitting lasso regression model using first alternative (iNim-RW) and MCMC

#Packages
library(nimble)
library(INLA)
library(mvtnorm)
library(MASS)
library(parallel)
library(coda)
library(ggmcmc)
library(myphdthesis)

# load data
load("bayesianLasso/hitters_data.RData")
df <- lassoDataDF

#Bayesian Lasso regression

fixedVals <- c("alpha","sigma")

#Define the conditional model to be fitted using INLA
fit.inla <- function(x ,
                     y ,
                     beta,
                     fixedVals,
                     #interInModel,
                     family
){

  data <- list(y=y, x=x)
  data$oset = data$x %*% matrix(beta, ncol = 1)
  res = INLA::inla(y ~ 1 + offset(oset),
                   data = data,
                   family = family,
                   control.fixed = list(prec.intercept = 0.001),
                   control.predictor = list(compute=TRUE))

  fitted_values = c(res$mlik[1,1])

  samples <- inla.posterior.sample(1, res)
  
  intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])


  ret <- data.frame(mld = fitted_values,
                    intercept,
                    precision,
                    row.names = NULL)
  colnames(ret) <- c("mld", fixedVals)

  ret <- as.matrix(ret)
  #ret <- c(ret)
  return(ret)
}

# compile NIMBLE function
nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1), # beta is a vector
    fixedVals = character(1, default = c("alpha","sigma")),
    #interInModel = double(0, default = 1),
    family = character(0, default = "gaussian")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)


# Define NIMBLE code
code <- nimbleCode({

  alpha ~ dnorm(0,0.01)

  for(j in 1:P){
    beta[j] ~ ddexp(location = 0, rate=est_lam)
  }

  #Fitting the inla with the simulated parameters
  for(i in 1:N){
  linpred[i] <- beta[1]*x[i,1]+beta[2]*x[i,2]+beta[3]*x[i,3]+beta[4]*x[i,4]+beta[5]*x[i,5] + alpha
  }

  sigma ~ dgamma(1,0.00005)

  # linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i],tau=sigma )
  }

})

## Parameterising the nimble model

#Data
df = lassoDataDF
inla_data <- list(y=as.numeric(df$y),
                  x = df$x,
                  y_obs=as.numeric(df$y))

#Constants
const <- list(N = length(df$y),
              P= ncol(df$x),
              est_lam = 1/0.0337
)

# Initial values
idm_inits <- function(){list(beta=rep(1,const$P),
                             sigma = 1
)
}

initsList <- idm_inits()

data = df

stdev.samp <- .25 * solve(t(data$x)%*%data$x)

samplers <- c("RW_INLA_block", "AFSS_INLA_block", "RW_block")
inlaMCMCtype <- c("inlamcmc", "inlamcmc", "mcmc")

#bayesianLasso <- list()

for(i in seq_along(samplers)){
  bayesianLasso <- INLAWiNimDataGenerating(data = c("y"),
                               covariate = data$x,
                               code = code,
                               family = "gaussian",
                               modelData = inla_data,
                               modelConstants = const,
                               modelInits = idm_inits,
                               nimbleINLA = nimbleINLA,
                               inlaMCMC = inlaMCMCtype[i],
                               inlaMCsampler = samplers[i],
                               samplerControl = list(propCov = stdev.samp,
                                                     #interInModel = 0,
                                                     #mu = c(rep(0,5)),
                                                    # scale = 1,
                                                     adaptive = FALSE),
                               parametersToMonitor = list(inla = c("alpha","sigma"),
                                                          mcmc = c("beta")),
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
  save(bayesianLasso, file = paste0("bayesianLasso/inlaWithNimble/inlaNimBayesianLasso",i,".RData"))
}




