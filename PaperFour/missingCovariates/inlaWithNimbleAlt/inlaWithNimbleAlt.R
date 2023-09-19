# Imputing missing covarites using second alternative (iNim2-RW) and MCMC

#load packages
library(mice)
data(nhanes2)
library(nimble)
library(inlamcmcnimble)

# data
d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi)) # finding na's
n.mis <- length(idx.mis) # number of nans
d.mis = data.frame(age = as.numeric(d.mis$age),
                   bmi = d.mis$bmi,
                   chl = d.mis$chl)%>%
  dplyr::arrange(bmi)

df = list(d.mis = d.mis, idx.mis = idx.mis)

ageModelMat <- model.matrix(~ as.factor(age),
                            data = d.mis,
                            contrasts.arg = list(age = contrasts(as.factor(d.mis$age), contrasts = TRUE)))

#model matrix covariates for NIMBLE
modelMat <- cbind(ageModelMat[, 2:3], d.mis$bmi)

# Define conditional model to be fitted using INLA
fitINLAMissingValues <- function(x ,
                                 y ,
                                 beta){
  missData <- data.frame(age= x[,1],
                         age1 = x[,2],
                         bmi = x[,3],
                         chl = y)
  #subset missing values
  idx.mis <- c(which(is.na(missData[,3])))

  missData[idx.mis,3] = beta

#Fit INLA model
  res = inla(chl ~ 1 + age + age1 + bmi,
             data = missData,
             family = "gaussian",
             verbose=FALSE,
             control.fixed = list(prec.intercept = 0.001),
             control.compute = list(config = TRUE),
             control.predictor = list(compute = TRUE),
             control.inla = list(cmin = 0))

  res <- inla.rerun(res)

  samples <- inla.posterior.sample(1, res)

  fitted_values <- samples[[1]]$latent[grepl("Predictor",rownames(samples[[1]]$latent) ),1]
  intercept = samples[[1]]$latent[grepl("Intercept",rownames(samples[[1]]$latent) ),1]
  precision <-  samples[[1]]$hyperpar[1]
  estbeta1 <- samples[[1]]$latent[grepl("age:1",rownames(samples[[1]]$latent) ),1]
  estbeta2 <- samples[[1]]$latent[grepl("age1:1",rownames(samples[[1]]$latent) ),1]
  estbeta3 <- samples[[1]]$latent[grepl("bmi:1",rownames(samples[[1]]$latent) ),1]
  ret <- cbind(fitted_values,
               intercept,
               estbeta1,
               estbeta2,
               estbeta3 ,
               precision)
  return(ret)
}


nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1)#, # beta is a vector
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fitINLAMissingValues'
)

nimbleINLA(x = modelMat,
           y = d.mis$chl,
           beta = rep(30, 9))


code <- nimbleCode({

beta0 <- inla.res[1, 2]
beta1 <- inla.res[1, 3]
beta2 <- inla.res[1, 4]
beta3 <- inla.res[1, 5]

  for(i in 1:n.idx){
    eta[i] ~ dnorm(muMiss, var = covMiss)
  }

  #Fitting the inla with the simulated parameters
   inla.res[1:N, 1:6] <- nimbleINLA(x[1:N,1:3], yObs[1:N], eta[1: n.idx])

  sigma <- inla.res[1, 6]

  #for(i in 1:N){
    linpred[1:N] <- inla.res[1:N, 1]#beta0 + beta1*x[i, 1] + beta2*x[i,2] + beta3*xobs[i]
  #}

  # linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i], tau=sigma)
  }

})

inla_data <- list(y = as.numeric(d.mis$chl),
                  yObs = as.numeric(d.mis$chl),
                  x = modelMat
)

#Constants
const <- list(N = nrow(inla_data$x),
              n.idx = length(df$idx.mis),
              muMiss = mean(inla_data$x[,3], na.rm = T),
              var = (4*var(inla_data$x[,3], na.rm = T))
)

# Initial values
idm_inits <- function(){list(eta = rep(mean(inla_data$x[,3], na.rm = T), const$n.idx)
)
}

# Fit the INLA within Nimble model
samplers <- c('RW_block')
nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
missingCovsAlt <- list()
for(i in seq_along(samplers)){
  missingCovsAlt = INLAWiNim(data = df,
                                    code = code,
                                    modelData = inla_data,
                                    modelConstants = const,
                                    modelInits = idm_inits,
                                    fam = "gaussian",
                                    mcmcSamplerChange = TRUE,
                                    parametersForSamplerChange = "eta",
                                    parametersToMonitor = c("beta0", "beta1", "beta2", "beta3", "sigma", "eta"),
                                    newSampler = samplers[i],
                                    newSamplerControl = list(scale = sqrt(10),
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

save(missingCovsAlt, file = "missingCovariates/inlaWithNimbleAlt/inlaNimMissingCovariatesAlt.RData")




