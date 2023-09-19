# Imputing missing covarites using first alternative (iNim-RW) and MCMC

# Load packages
library(mice)
data(nhanes2)
library(nimble)
library(inlamcmcnimble)
library(INLA)

# data
d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi)) # finding na's
n.mis <- length(idx.mis) # number of nans
d.mis = data.frame(age = as.numeric(d.mis$age),
              bmi = d.mis$bmi,
              chl = d.mis$chl)%>%
  dplyr::arrange(bmi)

df = list(d.mis = d.mis, idx.mis = idx.mis)

ageModelMat <- model.matrix(~ -1 + as.factor(age),
             data = data.frame(age = d.mis$age))

#model matrix covariates for NIMBLE
modelMat <- cbind(ageModelMat[, 2:3], d.mis$bmi)

# Fit conditional model using INLA
fitINLAMissingValues <- function(x ,
                                 y ,
                                 beta,
                                 fixedVals,
                                 #interInModel,
                                 family){
  missData <- data.frame(age1= x[,1],
                         age2 = x[, 2],
                         bmi = x[,3],
                        # bmi = x[,3],
                         chl = y)
  #subset missing values

  #missData <- d.mis
  idx.mis <- c(which(is.na(missData[,3])))

  missData[idx.mis,3] = beta

  # prior.fixed <- list(initial = 1, prior = "loggamma",
  #                     param = c(1,0.00005), fixed = FALSE)

  res = inla(chl ~ 1 + age1 + age2 +bmi,
             data = missData,
             family = family,
             verbose=FALSE,
             control.fixed = list(prec.intercept = 0.001),
             control.predictor = list(compute = TRUE),
             control.inla = list(cmin = 0)
  )

  beta0 = res$marginals.fixed[[1]]
  beta1 = res$marginals.fixed[[2]]
  beta2 = res$marginals.fixed[[3]]
  tau = res$marginals.hyperpar[[1]]
  samples <- inla.posterior.sample(1, res)
  
  fitted_values <- samples[[1]]$latent[grepl("Predictor",rownames(samples[[1]]$latent) ),1]
  intercept = samples[[1]]$latent[grepl("Intercept",rownames(samples[[1]]$latent) ),1]
  precision <-  samples[[1]]$hyperpar[1]
  estbeta1 <- samples[[1]]$latent[grepl("age:1",rownames(samples[[1]]$latent) ),1]
  estbeta2 <- samples[[1]]$latent[grepl("age1:1",rownames(samples[[1]]$latent) ),1]
  estbeta3 <- samples[[1]]$latent[grepl("bmi:1",rownames(samples[[1]]$latent) ),1]
  
  # precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])
  # intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  # estbeta1 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[2]])
  # estbeta2 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[3]])
  # estbeta3 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[4]])
  fitted_values = c(res$mlik[1,1])
  #fitted_values = res$summary.fitted.values[,"mean"]
  ret <- cbind(fitted_values,
               intercept,
               estbeta1,
               estbeta2,
               estbeta3,
               precision)
  colnames(ret) <- c("mld", fixedVals)
  return(ret)
}

# Compile INLA function
nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1), # beta is a vector
    fixedVals = character(1, default = c("alpha", "beta[1]","beta[2]","beta[3]", "sigma")),
    #interInModel = double(0, default = 1),
    family = character(0, default = "gaussian")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fitINLAMissingValues'
)

# Define NIMBLE code
code <- nimbleCode({
 alpha ~ dnorm(0, tau = 0.1) #dunif(-1000, 1000)
  for( i in 1:3){
    beta[i] ~ dnorm(0, tau = 0.001)
  }
  #alpha ~ dnorm(0, 0.001)
for(i in 1:n.idx){
  eta[i] ~ dnorm(muMiss, var = covMiss)
}
  #eta[1: n.idx] ~ dmnorm(muMiss[1:n.idx], cov = covMiss[1:n.idx, 1:n.idx])

   for(i in 1:16){
     xobs[i] <- x[i,3]
   }
  #
   for(i in 17:25){
     xobs[i] <- eta[i - 16]
   }
  #Fitting the inla with the simulated parameters
 # inla.res[1:N, 1:5] <- nimbleINLAMissingValues(x[1:N,1:3], idxMiss[1:n.idx], eta[1: n.idx])

  sigma ~ dgamma(1,0.00005)
  for(i in 1:25){
  linpred[i] <- beta[1]*x[i, 1] + beta[2]*x[i,2] + beta[3]*xobs[i]
  }

  #for(i in 17:25){
   # linpred[i] <- beta[1] + beta[2]* x[i,1] + beta[3]* eta[i - 16]  + beta[4]*x[i,3]
  #}
    #inprod(beta[1:P], x[i, 1:P])


  # linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i], tau=sigma)
  }

})

inla_data <- list(y = as.numeric(d.mis$chl),
                  x = modelMat
)

#covMiss <- 4 * ()
#Constants
const <- list(N = nrow(inla_data$x),
              n.idx = length(df$idx.mis),
              #muMiss = rep(mean(inla_data$x[,3], na.rm = T), length(df$idx.mis)),
              #covMiss = diag(mean(inla_data$x[,3], na.rm = T), length(df$idx.mis))
              muMiss = mean(inla_data$x[,3], na.rm = T),
              covMiss = (4*var(inla_data$x[,3], na.rm = T)),
              age = inla_data$x[,1]
)

# Initial values
idm_inits <- function(){list(eta = rep(mean(inla_data$x[,3], na.rm = T), const$n.idx),
                             beta = rep(1, 3),
                             alpha = 30
)
}

x = cbind(d.mis$age, d.mis$bmi); y = as.numeric(d.mis$chl)
nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)
samplers <- c("RW_INLA_block", "AFSS_INLA_block", "RW_block")
inlaMCMCtype <- c("inlamcmc", "inlamcmc", "mcmc")

#missingCovs <- list()

for(i in seq_along(samplers)){
  missingCovs <- INLAWiNimDataGenerating(data = c("y"),
                               covariate = x,
                               code = code,
                               family = "gaussian",
                               modelData = inla_data,
                               modelConstants = const,
                               modelInits = idm_inits,
                               nimbleINLA = nimbleINLA,
                               inlaMCMC = inlaMCMCtype[i],
                               inlaMCsampler = samplers[i],
                               samplerControl = list(#interInModel = 0,
                                         #            mu = c(rep(mean(inla_data$x[,3], na.rm = T),9)),
                                                     scale = sqrt(10),
                                                     adaptive = FALSE),
                               parametersToMonitor = list(inla = c("alpha", "beta", "sigma"),
                                                          mcmc = c("eta")),
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
  save(missingCovs, file = paste0("missingCovariates/inlaWithNimble/inlaNimMissingCovariates",i,".RData"))
}


