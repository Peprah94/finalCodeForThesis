# load packages
library(AHMbook)
RNGversion("3.5.3")
library(INLA)
library(inlabru)
library(sp)
library(sf)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)
library(ggplot2)
library(gridExtra)
library(readr)
#library(terra)
library(tidyr)
library(stringr)
data("BerneseOberland")
library(nimble)
library(myphdthesis)

ii = 0

load("spatialOccupancyModel/dataSimulated.RData")

N <- dataSimulated$nsurveys
indxOfSitesObs <- which(complete.cases(dataSimulated$yobs) == TRUE)
#dataSimulated$yobs[, ]

dataForModel <- data.frame(Longitude = dataSimulated$xcoord[indxOfSitesObs],
                           Latitude = dataSimulated$ycoord[indxOfSitesObs],
                           elevationS = dataSimulated$elevationS[indxOfSitesObs],
                           obsY = dataSimulated$y[indxOfSitesObs,1],
                           i_year = as.integer(factor(rep(1, each = length(dataSimulated$xcoord[indxOfSitesObs]))) ))%>%
  as.matrix()

x <- cbind(as.matrix(dataForModel[, c(1,2,3,5)]),
           rep(standardize(c(dataSimulated$forest[indxOfSitesObs])), 1),
           c(dataSimulated$wind[indxOfSitesObs,1]))%>%
  as.matrix()

x[ ,1] <- x[ ,1]/1000
x[ ,2] <- x[ ,2]/1000
x <- cbind(x, x[ ,3]^2)

# Define conditional model to be fitted using INLA
inlabruModelFit <- function(x, #matrix
                             y
                     ){

  # Simulate data
  coordsData <- x%>%
    as.data.frame()


  #convert y to vector
  y <- c(y)
  colnames(coordsData) <- c("Longitude","Latitude","elev","i_year" , "forest", "wind", "elevsq" )

  coordsDataNew <- sf::st_as_sf(coordsData,
                                coords = c("Longitude",
                                           "Latitude"))

  max.edge = diff(range(sf::st_coordinates(coordsDataNew)[,1]))/(3*8)
  bound <- diff(range(sf::st_coordinates(coordsDataNew)[,1]))/(3*5)

  mesh1 = inla.mesh.create.helper(points = st_coordinates(coordsDataNew),
                                  offset=c(1, 5),
                                  max.edge = c(10, 20),
                                  cutoff = 0)

  #SPDE with exponential correlation function
  spde = inla.spde2.matern(mesh = mesh1,
                           alpha = 1.5)

  cmp <- as.formula("obsY~ - 1 + beta0(1) + elev + elevsq + w2(main = coordinates, model = spde)")

  coordsDataNew1 <- as.data.frame(coordsData)

  coordsDataNew1 <- sp::SpatialPointsDataFrame(coords = coordsDataNew1[, c("Longitude", "Latitude")],
                                               data = coordsDataNew1)

  coordsDataNew1$obsY <- y

  lik1 <- inlabru::like("binomial",
                        formula = as.formula(paste0("obsY ~ beta0 + elev + elevsq +w2")),
                        data = coordsDataNew1,
                        Ntrials = 1
  )



  errorIndicator <- inherits(try( m_bru <- inlabru::bru(cmp,
                                                        lik1,
                                                        options =
                                                          list(
                                                            bru_verbose = TRUE,
                                                            bru_max_iter=1,
                                                            control.fixed = list(expand.factor.strategy = "inla",
                                                                                 mean = 0,
                                                                                 prec.intercept = 0.001,
                                                                                 prec = 1 / (100 * 100)),
                                                            control.family = list(link = "logit"),
                                                            control.predictor = list(link = 1),
                                                            control.inla = list(int.strategy = "eb",
                                                                                control.vb=list(enable=FALSE),
                                                                                cmin = 0.01),#, strategy = "gaussian"),
                                                            control.compute=list(return.marginals.predictor=TRUE)
                                                          )
  ), silent = TRUE),
  "try-error")
  #m_bru <- bru_rerun(m_bru)
  #indx <- mesh1$idx$loc
  # ret <- matrix(m_bru$summary.fitted.values[indx,"mean"], nrow = length(coordsData$Longitude)/N, ncol = N, byrow = FALSE)
  if(errorIndicator){
    mld = -Inf
    fittedValues <- -Inf
    intercept = NA
    elev = NA
    elevsq = NA
    siteSD = NA
    #siteSD <-  samples[[1]]$hyperpar[1]
    theta1 <-  NA
    theta2 <- NA
  }else{
    samples <- inla.posterior.sample(1, m_bru)
    fittedValues = c(m_bru$mlik[1,1])

    fitted_values <- samples[[1]]$latent[grepl("APredictor",rownames(samples[[1]]$latent) ),1]
    fittedVals <- fitted_values

    intercept = samples[[1]]$latent[grepl("beta0:1",rownames(samples[[1]]$latent) ),1]
    elev = samples[[1]]$latent[grepl("elev:1",rownames(samples[[1]]$latent) ),1]
    elevsq = samples[[1]]$latent[grepl("elevsq:1",rownames(samples[[1]]$latent) ),1]

    #siteSD <-  samples[[1]]$hyperpar[1]
    theta1 <-  samples[[1]]$hyperpar[1]
    theta2 <- samples[[1]]$hyperpar[2]
    # intercept, elev,
    ret1 <- cbind(fittedValues,intercept, elev, elevsq, theta1,  theta2,fittedVals)
  }
  #ret <- as.matrix(ret)
  return(ret1)

}

nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1)#, #y is a matrix
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'inlabruModelFit'
)

#nimbleINLA(x, dataSimulated$y[indxOfSitesObs, 1])

# Define distributions
dSpatOcc <- nimbleFunction(
  run = function(x = double(1), initial = double(1), covariate = double(2),
                 log = logical(0, default = 0)) {
    returnType(double())
    #return_marginal likelihood
runINLA <- nimbleINLA(x = covariate, y = x)
  mld <- runINLA[1,1]
    if(log) return(mld)
    return(exp(mld))
  })

rSpatOcc <- nimbleFunction(
  run = function(n = integer(), initial = double(1), covariate = double(2)) {
    returnType(double(1))
    isNA <- initial == 0
    nVars <- length(initial)
runINLA <- nimbleINLA(x = covariate, y = initial)
   fittedVals <- plogis(runINLA[1:nVars,7])
   samplesZ <- rbinom(nVars, 1, fittedVals)
    initial[isNA] <- samplesZ[isNA]

    #print(initial)
    return(initial)
  })

registerDistributions(list(
  dSpatOcc = list(
    BUGSdist = "dSpatOcc(initial, covariate)",
    discrete = TRUE,
    mixedSizes = TRUE,
    #range    = c('lower = 0', 'upper = 1')#,
    range = c(0, 1),
    types = c('value = double(1)', 'initial = double(1)', 'covariate = double(2)')
  )))

# NIMBLE code
code <-nimbleCode({
  # Specify priors
  for(v in 1:3){
    alpha[v] ~ dnorm(0, 0.001)
  }

  for(site.tag in 1:nsites){
    for(visit.tag in 1:nvisits){
      logit(p[site.tag, visit.tag]) <- alpha[1] + alpha[2]*forest[site.tag] + alpha[3]*wind[site.tag, visit.tag]
    }
  }

  beta0 <- inlaFit[1, 2]
  beta1 <- inlaFit[1, 3]
  beta2 <- inlaFit[1, 4]
  #siteSD <- inlaFit[1, 5]
  theta1 <- inlaFit[1, 5]
  theta2 <- inlaFit[1, 6]


  ###### Occupancy model

 #logit(psi[1:nsites]) <- inlaFit[1:nsites, 1]

   inlaFit[1:nsites, 1: 7] <- nimbleINLA(xObs[1:50, 1:7], z[1:nsites])#, alpha[1:3])

  # Observation model

   #for(site.tag in 1:nsites){
     z[1:nsites] ~ dSpatOcc(initial = zObs[1:nsites], covariate = xObs[1:50, 1:7])
   #}

   # Observation model
   for(site.tag in 1:nsites){
     for(visit.tag in 1:nvisits){
       y[site.tag, visit.tag] ~ dbin(size = 1, prob = z[site.tag]*p[site.tag, visit.tag])
     }
   }

})

zst <- apply(dataSimulated$y[indxOfSitesObs,], 1, function(t){
  r <-sum(t)
  rt <- as.numeric(r>0)
  #rt <- ifelse(rt==0, NA, rt)
  return(rt)
})

#nimbleINLA(x, zst, c(1, -1, 2))

## Parameterising the nimble model

#Data
inla_data <- list(y = dataSimulated$y[indxOfSitesObs, ],
                  yObs = dataSimulated$y[indxOfSitesObs, ],
                  forest = standardize(dataSimulated$forest[indxOfSitesObs]),
                  wind = dataSimulated$wind[indxOfSitesObs,],
                  elev = dataSimulated$elevationS[indxOfSitesObs],
                  xObs = x ,
                  elevsq = (dataSimulated$elevationS[indxOfSitesObs])^2,
                  zObs = zst)

#Constants
const <- list(N = dataSimulated$nsurveys,
              nvisits = dataSimulated$nsurveys,
              nsites = length(dataSimulated$xcoord[indxOfSitesObs])
)

# Initial values
idm_inits <- function(){list(alpha = c(1,-1, 1),
                             p = matrix(runif(const$N * const$nsites, 0, 1), nrow = const$nsites, const$N),
                             psi = runif(const$nsites, 0, 1),
                             z = zst
)
}

initsList <- idm_inits()


# Fit the INLA within Nimble model
samplers <- c("AF_slice_binary")
spatialModelAlt <- list()
for(i in seq_along(samplers)){
  spatialModelAlt = INLAWiNim(data = c("z"),
                                       code = code,
                                       modelData = inla_data,
                                       modelConstants = const,
                                       modelInits = idm_inits,
                                       fam = "binomial",
                                       mcmcSamplerChange = TRUE,
                                       parametersForSamplerChange = c("z"),
                                       parametersToMonitor = c("beta0",
                                                               "beta1",
                                                               "beta2",
                                                               #"siteSD",
                                                               "theta1",
                                                               "theta2",
                                                               "alpha",
                                                               "z"),
                                       newSampler = samplers[i],#'AF_slice',
                                       newSamplerControl = list(maxContractions = 3,
                                                                sliceMaxSteps = 3),
                                       mcmcConfiguration =  list(n.chains = 1,
                                                                 n.iterations = 1000,#0500,#500, #5500,#500,
                                                                 n.burnin = 200,#0500,#500,#500,#500,
                                                                 n.thin = 1,
                                                                 setSeed = TRUE,
                                                                 samples=TRUE,
                                                                 samplesAsCodaMCMC = TRUE,
                                                                 summary = TRUE,
                                                                 WAIC = FALSE))
}

save(spatialModelAlt, file = "spatialOccupancyModel/inlaWithNimbleAlt/spatialOccupancyAltNew.RData")

