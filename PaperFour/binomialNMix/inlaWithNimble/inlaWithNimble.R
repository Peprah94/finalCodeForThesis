# Fitting binomial N-Mixture model using first alternative (iNim-RW) and MCMC

library(unmarked)
data("mallard")
library(myphdthesis)
library(nimble)
library(INLA)
library(inlabru)

# Extract covariates
length <- mallard.site[ , "length"]
elev <- mallard.site[, "elev"]
forest <- mallard.site[, "forest"]
mean.ivel <- mallard.obs$ivel
mean.date <- mallard.obs$date

# Fit N-Mixture model with unmarked package
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

# Define confitional model to be fitted using INLA
fit.inla <- function(x ,
                     y ,
                     beta,
                     extraVars,
                     fixedVals,
                     family
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
  # samples <- inla.posterior.sample(1, res)
  # fitted_values <- samples[[1]]$latent[grepl("Predictor",rownames(samples[[1]]$latent) ),1]
  # intercept = samples[[1]]$latent[grepl("(Intercept)",rownames(samples[[1]]$latent) ),1]
  # length = samples[[1]]$latent[grepl("length",rownames(samples[[1]]$latent) ),1]
  # elev = samples[[1]]$latent[grepl("elev",rownames(samples[[1]]$latent) ),1]
  # forest = samples[[1]]$latent[grepl("forest",rownames(samples[[1]]$latent) ),1]
 if(errorIndicated){
   ret <- data.frame(mld = -Inf,
                     NA,
                     NA,
                     NA,
                     NA,
                     # rho,
                     row.names = NULL)
 }else{
 fitted_values = c(res$mlik[1,1])
 samples <- inla.posterior.sample(1, res)
 intercept = samples[[1]]$latent[grepl("(Intercept)",rownames(samples[[1]]$latent) ),1]
 beta1 = samples[[1]]$latent[grepl("ivel",rownames(samples[[1]]$latent) ),1]
beta2 = samples[[1]]$latent[grepl("date",rownames(samples[[1]]$latent) ),1][1]
beta3 = samples[[1]]$latent[grepl("date.sq",rownames(samples[[1]]$latent) ),1]
 # intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
 # beta1 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[2]])
 # beta2 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[3]])
 # beta3 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[4]])
  # rho = INLA::inla.emarginal(function(x) x,res$marginals.hyperpar[[1]])
  #precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])


  ret <- data.frame(mld = fitted_values,
                    intercept,
                    beta1,
                    beta2,
                    beta3,
                   # rho,
                    row.names = NULL)
}
 # colnames(ret) <- c("mld", fixedVals)

  ret <- as.matrix(ret)
  #ret <- c(ret)
  return(ret)
}

# Compile the INLA function
nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(2), #y is a matrix
    beta=double(1), # beta is a vector
    extraVars = double(1),
    fixedVals = character(1, default = c("intercept", "beta1", "beta2", "beta3")),
    #interInModel = double(0, default = 1),
   family = character(0, default = "binomial")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)

#nimbleINLA(x=x, extraVars = c(N), y = y, beta = c(0,0,0, 0))

code <- nimbleCode({

  for(i in 1:nsites){
    log(lambda[i]) <- beta[1] + beta[2]*x[i,1] + beta[3]*x[i,2] + beta[4]*x[i,3]
  }
  #lambda[1:N] <- nimbleINLA(x[1:N, 1:9], Nobs[1:N])

  for(i in 1:nsites){
    N[i] ~ dpois(lambda[i])
  }

  for(i in 1:nsites){
    for(j in 1:nvisits){
      logit(p.tag[i,j]) <- alpha[1] + alpha[2]*x[i, 3 + j] + alpha[3]*x[i, 6 + j] + alpha[4]*pow(x[i, 6 + j],2)
    }
  }

   # for(i in 1:nsites){
   #  for(j in 1:nvisits){
   #     pNb[i] <- rho/ (rho + (lambda[i]))
   #  }
   # }


 # rho ~ dunif(0.0001, 10)

  for(i in 1:nsites){
    for(j in 1:nvisits){
    y[i,j] ~ dbin(prob = p.tag[i,j], size = N[i])
    }
  }

  #Prior distributions
  for(i in 1:4){
    alpha[i] ~ dnorm(0, 0.001)
  }

  for(i in 1:4){
    beta[i] ~ dnorm(0, 0.001)
  }


  # Derived quantity
  Ntotal <- sum(N[1:nsites])
})


inla_data <- list(y=mallard.inla.df[,1:3],
                  x = x,
                  y_obs=mallard.inla.df[,1:3])

#Constants
const <- list(nsites = nrow(inla_data$y),
              nvisits = ncol(inla_data$y)

)

# Initial values
idm_inits <- function(){list(beta = rep(1, 4),
                             alpha = rep(0, 4),
                             rho = 0.8,
                             N = rep(30, const$nsites)
)
}

initsList <- idm_inits()

samplers <- c("RW_INLA_blockV2","AFSS_INLA_blockV2",  "AF_slice")
inlaMCMCtype <- c("inlamcmc", "inlamcmc", "mcmc")

#bayesianNMix <- list()
#nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
for(i in seq_along(samplers)){
#i <- 1
  bayesianNMix <- INLAWiNimDataGeneratingTargetDivide(data = c("y"),
                                                covariate = as.matrix(x),
                                                code = code,
                                                family = "binomial",
                                                modelData = inla_data,
                                                modelConstants = const,
                                                modelInits = idm_inits,
                                                nimbleINLA = nimbleINLA,
                                                inlaMCMC = inlaMCMCtype[i],
                                                inlaMCsampler = samplers[i],
                                                samplerControl = list(#propCov = stdev.samp,
                                                                      #interInModel = 0,
                                                                      #mu = c(rep(0,5)),
                                                                      scale = 1,
                                                                      adaptive = FALSE,
                                                                      extraVars = c("N"),
                                                                      adaptiveInterval = 100),
                                                parametersToMonitor = list(inla = c("alpha"),
                                                                           mcmc = c("beta", "N"),
                                                                           mcmcINLA = c("alpha"),
                                                                           additionalPars = c("Ntotal")),
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
  save(bayesianNMix, file = paste0("binomialNMix/inlaWithNimble/binomialNMix", i, ".RData"))
}




