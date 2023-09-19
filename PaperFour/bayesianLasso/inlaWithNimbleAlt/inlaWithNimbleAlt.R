# load data
load("bayesianLasso/hitters_data.RData")
df <- lassoDataDF

#Packages
library(nimble)
library(INLA)
library(mvtnorm)
library(MASS)
library(parallel)
library(coda)
library(ggmcmc)
library(myphdthesis)

#ii <- 0
# nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
# nimbleOptions(MCMCsaveHistory = TRUE)

fit.inla <- function(x, y, beta){
 # ii  <- get("ii",envir =  parent.frame())
  #ii <- assign("ii",ii+1,envir = parent.frame())
  #print(ii)
  data <- list(y=y, x=x)
  data$oset = data$x %*% beta
  res = INLA::inla(y ~ 1 + offset(oset),
                   data = data,
                   control.compute = list(config = TRUE),
                   control.predictor = list(compute = TRUE))
  #res = INLA::inla.rerun(res)

  #generate samples
  samples <- inla.posterior.sample(1, res)
  fitted_values <- samples[[1]]$latent[grepl("Predictor",rownames(samples[[1]]$latent) ),1]
  intercept = samples[[1]]$latent[grepl("Intercept",rownames(samples[[1]]$latent) ),1]
  precision <-  samples[[1]]$hyperpar
  ret <- cbind(precision, intercept, fitted_values)
  return(ret)
}


nimbleINLA <- nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1) # beta is a vector
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)

# NIMBLE code
code <- nimbleCode({
  #Prior for beta1 and beta2
  alpha <- inla.res[1,2]

  for(j in 1:P){
    beta[j] ~ ddexp(location = 0, rate=est_lam)
  }

  #Fitting the inla with the simulated parameters
  inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:P],y_obs[1:N],beta[1:P])

  #linpred[1:N] <- inla.res[1:100,3]
  sigma <- inla.res[1,1]
  linpred[1:N] <-  inla.res[1:N, 3]

  # linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i],tau=sigma )
  }
})

## Parameterising the nimble model
data = df

stdev.samp <- .25 * solve(t(data$x)%*%data$x)

#Data
inla_data <- list(y=as.numeric(df$y),
                  x = df$x,
                  y_obs=as.numeric(df$y))

#Constants
const <- list(N = length(df$y),
              P= ncol(df$x),
              est_lam = 1/0.0337
              )

# Initial values
inits <- function(){list(beta=rep(0,const$P)
)
}

#initsList <- idm_inits()


# Fit the INLA within Nimble model
samplers <- c('RW_block')
bayesianLassoAlt <- list()
for(i in seq_along(samplers)){
  bayesianLassoAlt[[i]] = INLAWiNim(data = df,
                                       code = code,
                                       modelData = inla_data,
                                       modelConstants = const,
                                       modelInits = inits,
                                       fam = "gaussian",
                                       mcmcSamplerChange = TRUE,
                                       parametersForSamplerChange = "beta",
                                    parametersToMonitor = c("beta", "sigma", "alpha"),
                                       newSampler = samplers[i],#'AF_slice',
                                       newSamplerControl = list(propCov = stdev.samp,
                                                                #mu = c(rep(0,5)),
                                                                scale = 1,
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

save(bayesianLassoAlt, file = "bayesianLasso/inlaWithNimbleAlt/inlaNimBayesianLassoAlt.RData")

