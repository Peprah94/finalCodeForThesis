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


#simulate data
# dataSimulated <- simOccSpatial(nsurveys =  10,
#                                mean.psi = 0.9,
#                                beta = c(2 ,-2),
#                                mean.p = 0.4,
#                                alpha = c(-1, 1),
#                                sample.size =50,
#                                variance.RF = 1,
#                                theta.RF = 10,
#                                seeds = c(10, 100),
#                                show.plots = TRUE)
# save(dataSimulated, file = "spatialOccupancyModel/dataSimulated.RData")

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

#x0 <- seq(-1.3, 4.3, length = 50)
#y0 <- seq(-1.3,4.3, length = 50)
#gridlocs <- expand.grid(x0,y0)

x <- cbind(as.matrix(dataForModel[, c(1,2,3,5)]),
           rep(standardize(c(dataSimulated$forest[indxOfSitesObs])), 1),
           c(dataSimulated$wind[indxOfSitesObs,1]))%>%
  as.matrix()

x[ ,1] <- x[ ,1]/1000
x[ ,2] <- x[ ,2]/1000
x <- cbind(x, x[ ,3]^2)#coordsData$elevationS^2

inlabruModelFit1 <- function(x, #matrix
                            y, #matrix
                            beta, #parameters fitted from elsewhere, and should be a vector
                            fixedVals,
                            family){

  # Simulate data
  coordsData <- x%>%
    as.data.frame()


  #convert y to vector
  y <- c(y)

  #p <- plogis(beta[1] + beta[2]* x[,5] + beta[3]* x[,6])
  #p = c(beta)
  #coordsData <- cbind(coordsData, p)
  colnames(coordsData) <- c("Longitude","Latitude","elevationS","i_year" , "forest", "wind" )

  coordsDataNew <- sf::st_as_sf(coordsData,
                                  coords = c("Longitude",
                                                "Latitude"))

  max.edge = diff(range(sf::st_coordinates(coordsDataNew)[,1]))/(3*5)

  mesh1 = inla.mesh.2d(loc = st_coordinates(coordsDataNew),
                       max.edge = c(max.edge-10, max.edge+10))

  #SPDE with exponential correlation function
  spde = inla.spde2.matern(mesh = mesh1,
                             alpha = 1.5)

  elev.spix <- SpatialPixelsDataFrame(point =  st_coordinates(coordsDataNew),
                                      data = data.frame(elev = coordsDataNew$elevationS))

  elevsq.spix <- SpatialPixelsDataFrame(point =  st_coordinates(coordsDataNew),
                                        data = data.frame(elev = (coordsDataNew$elevationS)^2))

  cmp <- as.formula("obsY~ - 1 + beta0(1) + elev(main = elev.spix, model = 'linear') + elevsq(main = elevsq.spix, model = 'linear')+ site(main= i_year, model = 'iid', n = N) + w2(main = coordinates, model = spde)")

  coordsDataNew1 <- as.data.frame(coordsData)

  coordsDataNew1 <- sp::SpatialPointsDataFrame(coords = coordsDataNew1[, c("Longitude", "Latitude")],
                                               data = coordsDataNew1)

  coordsDataNew1$obsY <- y

  functionConstants <- function(p, intercept, elev, elevsq, sites, w2){
    #linearPred <- plogis(intercept + elev + elevsq + sites + w)
   linearPred <- intercept + elev + elevsq + sites + w2
   p <- abs(p-0.0001)
   #firstTerm <- log(p)/(1- log(1 + exp(linearPred)))
   #secondTerm <- linearPred / (1- log(1 + exp(linearPred)))
   firstTerm <- exp(linearPred) / (1+ exp(linearPred))
   secondTerm <- -log(p)
   denom <- linearPred -log(p) - log(1 + exp(linearPred))
    #
    #detProb <- qlogis(p)
    #firstTerm <- log( 1- linearPred*p)
    #secondTerm <- log(1 - linearPred)
    #thirdTerm <- log(1 - p)
   # ret <- detProb - firstTerm + secondTerm + thirdTerm
   denom <- denom - 0.001
   ret <- denom/(1-denom) #1- (firstTerm + log(secondTerm))#firstTerm + secondTerm
   return(ret)
  }

  functionConstants1 <- function(p, intercept, elev, elevsq, sites, w2){
    #linearPred <- plogis(intercept + elev + elevsq + sites + w)
    linearPred <- intercept + elev + elevsq + sites + w2
    p <- abs(p-0.0001)
   # firstTerm <- log(p)/(1- log(1 + exp(linearPred)))
   # secondTerm <- linearPred / (1- log(1 + exp(linearPred)))
    #firstTerm <- exp(linearPred) / (1+ exp(linearPred))
    #seconTerm <- -log(p)
    # p <- abs(p-0.0001)
    #detProb <- qlogis(p)
    #firstTerm <- log( 1- linearPred*p)
    #secondTerm <- log(1 - linearPred)
    #thirdTerm <- log(1 - p)
    # ret <- detProb - firstTerm + secondTerm + thirdTerm
    ret <- -log(p) -log(1 + exp(linearPred))#1- (firstTerm + secondTerm)#firstTerm + secondTerm
    return(ret)
  }

  lik1 <- inlabru::like(family,
                        formula = as.formula(paste0("obsY ~ beta0 + elev + elevsq + site +w2 ")),
                        #formula = as.formula(paste0("obsY ~ functionConstants(p, beta0, elev, elevsq, site, w2)")),
                       # formula = as.formula(paste0("obsY ~   beta0, elev, elevsq, site, w2)")),
                        data = coordsDataNew1,
                        Ntrials = 1,
                        domain = list(coordinates = mesh1)
  )



  m_bru <- inlabru::bru(cmp,
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
                            control.inla = list(int.strategy = "eb"),#, strategy = "gaussian"),
                            control.compute=list(return.marginals.predictor=TRUE)
                          )
  )
    #m_bru <- bru_rerun(m_bru)
  #indx <- mesh1$idx$loc
 # ret <- matrix(m_bru$summary.fitted.values[indx,"mean"], nrow = length(coordsData$Longitude)/N, ncol = N, byrow = FALSE)
  samples <- inla.posterior.sample(1, m_bru)
  fittedValues = c(m_bru$mlik[1,1])
  fitted_values <- samples[[1]]$latent[grepl("APredictor",rownames(samples[[1]]$latent) ),1]

  # fittedVals <- matrix(fitted_values,
  #                      nrow = length(coordsData$Longitude)/N,
  #                      ncol = N,
  #                      byrow = FALSE)
  fittedVals <- fitted_values

  intercept = samples[[1]]$latent[grepl("beta0:1",rownames(samples[[1]]$latent) ),1]
  elev = samples[[1]]$latent[grepl("elev:1",rownames(samples[[1]]$latent) ),1]
  elevsq = samples[[1]]$latent[grepl("elevsq:1",rownames(samples[[1]]$latent) ),1]

  siteSD <-  samples[[1]]$hyperpar[1]
  theta1 <-  samples[[1]]$hyperpar[2]
  theta2 <- samples[[1]]$hyperpar[3]

  ret1 <- cbind(fittedValues, intercept, elev, elevsq, theta1, siteSD, theta2)
  colnames(ret1) <- c("mld", fixedVals, "siteSD", "theta1")

  #ret <- as.matrix(ret)
  return(ret1)

}

inlabruModelFit <- function(x, #matrix
                             y, #matrix
                             beta, #parameters fitted from elsewhere, and should be a vector
                             fixedVals,
                             family){

  # Simulate data
  coordsData <- x%>%
    as.data.frame()


  #convert y to vector
  y <- c(y)

  #p <- plogis(beta[1] + beta[2]* x[,5] + beta[3]* x[,6])
  #p = c(beta)
  #coordsData <- cbind(coordsData, p)
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

  #elev.spix <- SpatialPixelsDataFrame(point =  st_coordinates(coordsDataNew),
   #                                   data = data.frame(elev = coordsDataNew$elevationS))

 # elevsq.spix <- SpatialPixelsDataFrame(point =  st_coordinates(coordsDataNew),
   #                                     data = data.frame(elev = (coordsDataNew$elevationS)^2))

  cmp <- as.formula("obsY~ - 1 + beta0(1) + elev + elevsq + w2(main = coordinates, model = spde)")

  coordsDataNew1 <- as.data.frame(coordsData)

  coordsDataNew1 <- sp::SpatialPointsDataFrame(coords = coordsDataNew1[, c("Longitude", "Latitude")],
                                               data = coordsDataNew1)

  coordsDataNew1$obsY <- y

  lik1 <- inlabru::like(family,
                        formula = as.formula(paste0("obsY ~ beta0 + elev + elevsq +w2")),
                        #formula = as.formula(paste0("obsY ~ functionConstants(p, beta0, elev, elevsq, site, w2)")),
                        #formula = as.formula(paste0("obsY ~  functionConstants(p, beta0, elev, elevsq, site, w2)")),
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
  # intercept = INLA::inla.emarginal(function(x) x,m_bru$marginals.fixed[[1]])
  # elev = INLA::inla.emarginal(function(x) x,m_bru$marginals.fixed[[2]])
  # elevsq = INLA::inla.emarginal(function(x) x,m_bru$marginals.fixed[[3]])
  # siteSD = INLA::inla.emarginal(function(x) x,1/m_bru$marginals.hyper[[1]])
  # theta1 = INLA::inla.emarginal(function(x) x,1/m_bru$marginals.hyper[[1]])
  # theta2 = INLA::inla.emarginal(function(x) x,1/m_bru$marginals.hyper[[2]])

  fitted_values <- samples[[1]]$latent[grepl("APredictor",rownames(samples[[1]]$latent) ),1]

  # fittedVals <- matrix(fitted_values,
  #                      nrow = length(coordsData$Longitude)/N,
  #                      ncol = N,
  #                      byrow = FALSE)
  fittedVals <- fitted_values

  intercept = samples[[1]]$latent[grepl("beta0:1",rownames(samples[[1]]$latent) ),1]
  elev = samples[[1]]$latent[grepl("elev:1",rownames(samples[[1]]$latent) ),1]
  elevsq = samples[[1]]$latent[grepl("elevsq:1",rownames(samples[[1]]$latent) ),1]

  #siteSD <-  samples[[1]]$hyperpar[1]
  theta1 <-  samples[[1]]$hyperpar[1]
  theta2 <- samples[[1]]$hyperpar[2]
# intercept, elev,
  ret1 <- cbind(fittedValues,intercept, elev, elevsq, theta1,  theta2)
  colnames(ret1) <- c("mld", fixedVals,  "theta1")
}
  #ret <- as.matrix(ret)
  return(ret1)

}


nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a matrix
    beta=double(1), # beta is a vector
    fixedVals = character(1, default = c("intercept", "beta1","beta2", "sigma")),
    #interInModel = double(0, default = 1),
    family = character(0, default = "binomial")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'inlabruModelFit'
)

zst <- apply(dataSimulated$y[indxOfSitesObs,], 1, function(t){
  r <-sum(t)
  rt <- as.numeric(r>0)
  #rt <- ifelse(rt==0, NA, rt)
  return(rt)
})

#text if it works
nimbleINLA(x, zst, c(1, -1, 2))


# NIMBLE code
code <-nimbleCode({
  # Specify priors
  alpha0 ~ dnorm(0, 0.1)
  beta0 ~ dnorm(0, 0.1)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.1)
    beta[v] ~ dnorm(0, 0.1)
  }

  for(site.tag in 1:nsites){
    for(visit.tag in 1:nvisits){
      logit(p[site.tag, visit.tag]) <- alpha0 + alpha[1]*forest[site.tag] + alpha[2]*wind[site.tag, visit.tag]
    }
  }


###### Occupancy model
  for(site.tag in 1:nsites){
      logit(psi[site.tag]) <- beta0 + beta[1]*elev[site.tag] + beta[2]*elevsq[site.tag] + eta[site.tag]
  }

  for(site.tag in 1:nsites){
z[site.tag] ~ dbin(size = 1, prob = psi[site.tag])
  }

  eta[1:nsites] ~ dcar_normal(adj[1:16], weights[1:16], num[1:50], tau = sigma)

  # Observation model
  for(site.tag in 1:nsites){
    for(visit.tag in 1:nvisits){
      y[site.tag, visit.tag] ~ dbin(size = 1, prob = z[site.tag]*p[site.tag, visit.tag])
    }
  }

  #
sigma ~ dgamma(1, 0.0001)

theta1 ~ dunif(-1000, 1000)

})


## Parameterising the nimble model
library(spdep)
coordgrid <- cbind(dataSimulated$xcoord[indxOfSitesObs], dataSimulated$ycoord[indxOfSitesObs])
neigh <- dnearneigh(coordgrid,
                    d1 = 0,
                    d2 = sqrt(2)*1000 + 1)
winnb <- nb2WB(neigh)
str(winnb)

zst <- apply(dataSimulated$y[indxOfSitesObs,], 1, function(t){
  r <-sum(t)
  rt <- as.numeric(r>0)
  rt <- ifelse(rt==0, NA, rt)
  return(rt)
})

#Data
inla_data <- list(y = dataSimulated$y[indxOfSitesObs,],
                  forest = standardize(dataSimulated$forest[indxOfSitesObs]),
                  wind = dataSimulated$wind[indxOfSitesObs,],
                  elev = dataSimulated$elevationS[indxOfSitesObs],
                  elevsq = (dataSimulated$elevationS[indxOfSitesObs])^2,
                  adj = winnb$adj,
                  weights = winnb$weights,
                  num = winnb$num,
                  z = zst)

#Constants
const <- list(N = dataSimulated$nsurveys,
              nvisits = dataSimulated$nsurveys,
              nsites = length(dataSimulated$xcoord[indxOfSitesObs])
)
#zst <- apply(inla_data$y, 1, max)
#zst[is.na(zst)] <- 1
# Initial values
idm_inits <- function(){list(p = matrix(runif(const$N * const$nsites, 0, 1), nrow = const$nsites, ncol= const$N),
                             alpha = c(-0.2, 0.2),
                             alpha0 = -0.405,
                             beta = c(0.1, 0.1),
                             beta0 = 1,
                             sigma = 1,
                             eta = rep(0, const$nsites),
                             psi = runif(const$nsites, 0, 1)
)
}

initsList <- idm_inits()

samplers <- c("RW_INLA_block", "RW_block")
inlaMCMCtype <- c("inlamcmc", "mcmc")
#samplers <- c("AFSS_INLA_block_binary", "AF_slice_binary")
#inlaMCMCtype <- c( "inlamcmc",  "mcmc")

#occSpatialModel <- list()
for(i in seq_along(samplers)){
  occSpatialModel <- INLAWiNimDataGeneratingTargetDivide(data = c("z"),
                               covariate = x,
                               code = code,
                               family = "binomial",
                               modelData = inla_data,
                               modelConstants = const,
                               modelInits = idm_inits,
                               nimbleINLA = nimbleINLA,
                               inlaMCMC = inlaMCMCtype[i],
                               inlaMCsampler = samplers[i],
                               samplerControl = list(maxContractions = 500,
                                                     sliceMaxSteps = 10),
                               parametersToMonitor = list(mcmc = c("alpha0","alpha", "z"),
                                                          mcmcINLA = c("beta0", "beta[1]","beta[2]", "sigma"),
                                                          inla = c("beta0", "beta[1]","beta[2]", "sigma")),
                               mcmcConfiguration = list(n.chains = 1,
                                                        n.iterations = 1000,#500,#500,
                                                        n.burnin = 200,#500,#00,
                                                        n.thin = 1,
                                                        setSeed = TRUE,
                                                        samples=TRUE,
                                                        samplesAsCodaMCMC = TRUE,
                                                        summary = TRUE,
                                                        WAIC = FALSE)
)
  save(occSpatialModel, file = paste0("spatialOccupancyModel/inlaWithNimble/occSpatialModelSmallNew",i,".RData"))
}

