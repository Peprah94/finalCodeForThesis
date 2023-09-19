library(unmarked)
data("mallard")
library(myphdthesis)
library(nimble)
library(INLA)
library(inlabru)


length <- mallard.site[ , "length"]
elev <- mallard.site[, "elev"]
forest <- mallard.site[, "forest"]
mean.ivel <- mallard.obs$ivel
mean.date <- mallard.obs$date

# unmarked
mallard.umf <- unmarkedFramePCount(y = mallard.y,
                                   siteCovs = mallard.site,
                                   obsCovs = mallard.obs)

out.unmk.2 <- pcount(~ ivel + date + I(date^2) ~ length + elev + forest,
                     mixture = "NB",
                     data = mallard.umf)
summary(out.unmk.2)

mallard.inla.df <- data.frame(y1 = mallard.y[ , "y.1"],
                              y2 = mallard.y[ , "y.2"],
                              y3 = mallard.y[ , "y.3"],
                              length = length,
                              elev = elev,
                              forest = forest,
                              mean.ivel,
                              mean.date)

#na.omit(mallard.inla.df[ , -c(1,2,3)])
mallard.inla.df <- mallard.inla.df[complete.cases(mallard.inla.df),]
fixedVals = c("intercept", "beta1", "beta2", "beta3")

x <- mallard.inla.df[ , -c(1,2,3)]
y <- mallard.inla.df[ , c(1,2,3)]


fit.inla <- function(x ,
                     y,
                     extraVars
){


  #convert y to vector
  n <- nrow(y)
  y <- c(unlist(y))
  ivel <- c(unlist(x[,4:6]))
  date <- c(unlist(x[,7:9]))
  extraVars <- as.numeric(extraVars)
  date.sq <- date^2
  site <- c(rep(1, n), rep(2,n), rep(3,n))
  data <- data.frame(ivel = ivel,
                     date = date,
                     date.sq = date.sq,
                     site = site,
                     y = y,
                     nTrials = rep(extraVars, 3))

  #p <- plogis(beta[1] + beta[2]*ivel + beta[3]* date + beta[4]* (date)^2)
  #list(y=y, x=x)
  #data$oset = data$x %*% beta
  errorIndicated <-  inherits(try(res <- INLA::inla(y ~ 1 + ivel + date + date.sq + f(site, model = "iid"),
                                                    data = data,
                                                    family = "binomial",
                                                    Ntrials =  nTrials,
                                                    control.family= list(link = 'logit'),
                                                    control.compute = list(config = TRUE),
                                                    control.predictor = list(compute=TRUE,
                                                                             link = 1)),
                                  silent = TRUE),
                              "try-error"
  )
  #generate samples
 if(errorIndicated){
   fitted_values <- rep(0, n)
   intercept = NA
   length = NA
   elev = NA
   forest = NA
   mld = -Inf
 }else{
   samples <- inla.posterior.sample(1, res)
   fitted_values <- plogis(samples[[1]]$latent[grepl("Predictor",rownames(samples[[1]]$latent) ),1])
   intercept = samples[[1]]$latent[grepl("(Intercept)",rownames(samples[[1]]$latent) ),1]
   ivel = samples[[1]]$latent[grepl("ivel",rownames(samples[[1]]$latent) ),1]
   date = samples[[1]]$latent[grepl("date",rownames(samples[[1]]$latent) ),1][1]
   date.sq = samples[[1]]$latent[grepl("date.sq",rownames(samples[[1]]$latent) ),1]
   mld = c(res$mlik[1,1])
 }


  # intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  # beta1 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[2]])
  # beta2 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[3]])
  # beta3 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[4]])
  # rho = INLA::inla.emarginal(function(x) x,res$marginals.hyperpar[[1]])
  #precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])

  #
  # ret <- data.frame(mld = mld,
  #                   intercept,
  #                   length,
  #                   elev,
  #                   forest,
  #                   # rho,
  #                   row.names = NULL)
  ret <- cbind(mld,
               fitted_values,
               intercept,
               ivel,
               date,
               date.sq)

  # colnames(ret) <- c("mld", fixedVals)

  ret <- as.matrix(ret)
  #ret <- c(ret)
  return(ret)
}


nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(2),
    extraVars = double(1)
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)

# Define distributions
dNmix <- nimbleFunction(
  run = function(x = double(2), extraVars = double(1), covariate = double(2), initial = double(2),
                 log = logical(0, default = 0)) {
    returnType(double())

    #return_marginal likelihood
    runINLA <- nimbleINLA(x = covariate, y = x, extraVars = extraVars)
    mld <- runINLA[1,1]

    if(log) return(mld)
    return(exp(mld))
  })

inlaMCMCMatrix <- function(samples, nrow, ncol){
  ret <- matrix(samples, nrow, ncol, byrow = FALSE)
  return(ret)
}

nimbleinlaMCMCMatrix <- nimble::nimbleRcall(
  prototype = function(
    samples=double(1), #x is a matrix
    nrow=integer(0),
    ncol = integer(0)
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'inlaMCMCMatrix'
)

rNmix <- nimbleFunction(
  run = function(n = integer(), initial = double(2), covariate = double(2), extraVars = double(1)) {
    returnType(double(2))
    #isNA <- initial == 0
    nVars <- length(extraVars)*3
    runINLA <- nimbleINLA(x = covariate, y = initial, extraVars = extraVars)
    fittedVals <- runINLA[1:nVars,1]

    #for(i in 1)
    samplesZ <- rbinom(nVars, prob = fittedVals, size = rep(extraVars, 3))
    ret <- nimbleinlaMCMCMatrix(samplesZ[1:nVars], nrow = nVars, ncol = 3)
      #matrix(samplesZ[1:nVars], nrow = nVars, ncol = 3)
    #initial[isNA] <- samplesZ[isNA]
    print(ret)
    return(initial)
  })

registerDistributions(list(
  dNmix = list(
    BUGSdist = "dNmix(initial, covariate, extraVars)",
    discrete = TRUE,
    mixedSizes = TRUE,
    #range    = c('lower = 0', 'upper = 1')#,
    range = c(0, Inf),
    types = c('value = double(2)', 'initial = double(2)', 'covariate = double(2)', 'extraVars = double(1)')
  )))



code <- nimbleCode({

  # for(i in 1:nsites){
  #   log(lambda[i]) <- beta[1] + beta[2]*x[i,1] + beta[3]*x[i,2] + beta[4]*x[i,3]
  # }
  #lambda[1:N] <- nimbleINLA(x[1:N, 1:9], Nobs[1:N])

  # for(i in 1:nsites){
  #   N[i] ~ dpois(lambda[i])
  # }

  inlaFit[1:573, 1: 6] <- nimbleINLA(x = x[1:nsites, 1:9], y = y_obs[1:nsites,1:3],extraVars = N[1:nsites])
  alpha[1] <- inlaFit[1, 3]
  alpha[2] <- inlaFit[1, 4]
  alpha[3] <- inlaFit[1, 5]
  alpha[4] <- inlaFit[1, 6]

  for(i in 1:nsites){
    log(lambda[i]) <- beta[1] + beta[2]*x[i,1] + beta[3]*x[i,2] + beta[4]*x[i,3]
  }
  #lambda[1:N] <- nimbleINLA(x[1:N, 1:9], Nobs[1:N])

  for(i in 1:nsites){
    N[i] ~ dpois(lambda[i])
  }

  for(i in 1:4){
    beta[i] ~ dnorm(0, 0.001)
  }
  # for(i in 1:nsites){
  #  for(j in 1:nvisits){
  #     pNb[i] <- rho/ (rho + (lambda[i]))
  #  }
  # }


  # rho ~ dunif(0.0001, 10)

  #for(i in 1:nsites){
   # for(j in 1:nvisits){
      y[1:nsites,1:3] ~ dNmix(initial = y_obs[1:nsites,1:3], covariate = x[1:nsites, 1:9], extraVars = N[1:nsites])
   # }
  #}

  #Prior distributions
  # for(i in 1:4){
  #   alpha[i] ~ dnorm(0, 0.001)
  # }



  # Derived quantity
  Ntotal <- sum(N[1:nsites])
})

Nst <- apply(mallard.inla.df[,1:3], 1, max)
Nst <- Nst + 10

inla_data <- list(y=mallard.inla.df[,1:3],
                  x = x,
                  nObs = Nst,
                  y_obs=mallard.inla.df[,1:3])

#Constants
const <- list(nsites = nrow(inla_data$y),
              nvisits = ncol(inla_data$y)

)

# Initial values
idm_inits <- function(){list(alpha = rep(0, 4),
                             beta = rep(0,4),
                             rho = 0.8,
                             N = Nst
)
}

initsList <- idm_inits()

samplers <- c( "RW_block")
inlaMCMCtype <- c("inlamcmc")

bayesianNMixAlt <- list()
#nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
for(i in seq_along(samplers)){
  bayesianNMixAlt <- INLAWiNim(data = x,
                                 code = code,
                                 modelData = inla_data,
                                 modelConstants = const,
                                 modelInits = idm_inits,
                                 fam = "binomial",
                                 mcmcSamplerChange = TRUE,
                                 parametersForSamplerChange = c("beta"),
                                 parametersToMonitor = c("beta",
                                                         "alpha",
                                                         "N",
                                                         "Ntotal"),
                                 newSampler = samplers[i],#'AF_slice',
                                 newSamplerControl = list(adaptive = FALSE),
                                 mcmcConfiguration =  list(n.chains = 1,
                                                           n.iterations = 100500,#500,#500, #5500,#500,
                                                           n.burnin = 50500,#500,#500,#500,
                                                           n.thin = 5,
                                                           setSeed = TRUE,
                                                           samples=TRUE,
                                                           samplesAsCodaMCMC = TRUE,
                                                           summary = TRUE,
                                                           WAIC = FALSE))
}

save(bayesianNMixAlt, file = "binomialNMix/inlaWithNimbleAlt/binomialNMixAlt.RData")







