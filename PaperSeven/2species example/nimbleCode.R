load("simData.RData")
load("detData.RData")
load("allDataForModelNew.RData")
detData <- detdata
source("DataGeneration.R")
library(sf)
library(terra)
library(raster)
library(inlabru)

for(iter in 1: 100){
  library(sf)
  library(terra)
  library(raster)
  library(inlabru)
  nspecies <- 2
  BNGproj <- CRS("+proj=robin +datum=WGS84 +units=km")
  cov1.rast <- simulateddata[[iter]]$cov[[1]]
  crs(cov1.rast) <- BNGproj
  cov1.spix <- as(cov1.rast,"SpatialPixelsDataFrame")
  input <- simulateddata[[iter]]$input
  #Covariates for first thinning
  # cov2.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate_thin.im$v))))))
  # r <- raster(cov2.sp)
  # r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
  # cov2.rast <- rasterize(cov2.sp@coords,r1,cov2.sp$cov, fun=mean,na.rm=T)
  cov2.rast <- simulateddata[[iter]]$cov[[2]]
  crs(cov2.rast) <- BNGproj
  cov2.spix <- as(cov2.rast,"SpatialPixelsDataFrame")

  #Covariate for second thinning
  # cov3.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate_detect.im$v))))))
  # r <- raster(cov3.sp)
  # r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
  cov3.rast <- simulateddata[[iter]]$cov[[3]] #rasterize(cov3.sp@coords,r1,cov3.sp$cov, fun=mean,na.rm=T)
  crs(cov3.rast) <- BNGproj
  cov3.spix <- as(cov3.rast,"SpatialPixelsDataFrame")

  ## Extra information on species detection ##




  ## Fit the model using inlabru ##

  ## the borders of the study region
  #coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
  #poly <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))

  ## the mesh
  #mesh <- inla.mesh.2d(loc.domain = coordsmat, offset = c(0.3, 1),
  #         max.edge = c(0.1, 0.5), cutoff = 0.2)
  # mesh <- simulateddata[[iter]]$mesh
  # ## SPDEs definition
  # spdes <- list()
  # for(i in 1: nspecies){
  #   spdes[[i]] <- inla.spde2.pcmatern(mesh = mesh,
  #                                     # PC-prior on range: P(practic.range < 0.05) = 0.01
  #                                     prior.range = c(input$ecological$hyperparameters$range[i], 0.1),
  #                                     # PC-prior on sigma: P(sigma > 1) = 0.01
  #                                     prior.sigma = c(sqrt(input$ecological$hyperparameters$sigma2[i]), 0.1))
  # }
  # 
  # #SPDEs for the thinning
  # spde2 <- inla.spde2.pcmatern(mesh = mesh,
  #                              # PC-prior on range: P(practic.range < 0.05) = 0.01
  #                              prior.range = c(input$sampling$hyperparameters$range, 0.01),
  #                              # PC-prior on sigma: P(sigma > 1) = 0.01
  #                              prior.sigma = c(sqrt(input$sampling$hyperparameters$sigma2), 0.01))

  csdata = simulateddata[[iter]]$thirdstage
  cssampdata = simulateddata[[iter]]$firststage$Samp_PPFinal
  detdata = detData[[iter]]
  covslist <- list(cov1.spix,cov2.spix,cov3.spix)
  #spdeslist <- list(spdes=spdes,spde2=spde2)
  covs = covslist
  region=simulateddata[[iter]]$region
  mesh=simulateddata[[iter]]$mesh

  data_df <- data.frame(
    Y = csdata$classifications$error,
    C = csdata$classifications$true_species,
    eco_cov = raster::extract(cov1.rast,csdata$classifications),
    samp_cov= raster::extract(cov2.rast,csdata$classifications),
    det_cov = raster::extract(cov3.rast,csdata$classifications))



  tmp <- csdata$classification
  # Eco_PPFinal_detect <- list()
  # for(i in 1:nspecies){Eco_PPFinal_detect[[i]] <- tmp[which(tmp$error==i),]
  # # fm_crs(Eco_PPFinal_detect[[i]]) <- fm_crs(covs[[1]])
  # slot(Eco_PPFinal_detect[[i]], "proj4string") <-  BNGproj
  # slot(detdata[[i]], "proj4string") <-  BNGproj
  # }
  # slot(cssampdata, "proj4string") <-  BNGproj
  source("estpar.R")

nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    p11=double(0), #x is a matrix
    p22 = double(0),
    model =character(0, "VSEDetect"),
    return = character(0, "inla")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'est_par'
)

# Define distributions
dSpatial <- nimbleFunction(
  run = function(x = double(1), p11 = double(0), p22 = double(0),initial = double(1),
                 log = logical(0, default = 0)) {
    returnType(double())

    #return_marginal likelihood
    #n <- nrow(omega)
    runINLA <- nimbleINLA(p11, p22, model = "VSEDetect", return = "inla")
    mld <- runINLA[1,4]

    if(log) return(mld)
    return(exp(mld))
  })

rSpatial <- nimbleFunction(
  run = function(n = integer(), p11 = double(0), p22 = double(0), initial = double(1)) {
    returnType(double(1))
    #isNA <- initial == 0
    nVars <- length(initial)
    #n <- nrow(omega)
    runINLA <- nimbleINLA(p11, p22, model = "VSEDetect", return = "inla")
    samplesZ <- runINLA[1:nVars, 3]

    #initial[isNA] <- samplesZ[isNA]
    #print(initial)
    return(samplesZ)
  })

registerDistributions(list(
  dSpatial = list(
    BUGSdist = "dSpatial(p11, p22,initial)",
    discrete = TRUE,
    mixedSizes = TRUE,
    #range    = c('lower = 0', 'upper = 1')#,
    range = c(1, 2),
    types = c('value = double(1)', 'initial = double(1)', 'p11 = double(0)', 'p22 = double(0)')
  )))


library(nimble)
library(dirmult)

# Estimate probability P(true = i| reported = j)
propTrueClass <- function(p11, p22, response){
  omega <- matrix(c(p11, 1-p11, 1-p22, p22), 2,2, byrow = TRUE)
  runINLA <- nimbleINLA(p11, p22, model = "VSEDetect", return = "inla")
  trueProbs <- runINLA[,1:2]
  reportedProbs <- matrix(NA, nrow = nrow(trueProbs), ncol = 2)
  for(i in 1: nrow(trueProbs)){
    for(j in 1:2){
    reportedProbs[i, j] <- (omega[response[i],j]*trueProbs[i, j])
    }
  }
  reportedProbs <- proportions(reportedProbs, margin = 1)
  return(reportedProbs)
}

nimblePropTrueClass <- nimble::nimbleRcall(
  prototype = function(
    p11=double(0), #x is a matrix
    p22 = double(0),
    response = double(1)
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'propTrueClass'
)

code <-nimbleCode({
  #prior for omega
  #for(i in 1:nspecies){
  #  alpha[i] ~ dexp(1)
  # }

  #for(i in 1:nspecies){
   # omega[i, 1:nspecies] ~ ddirch(alpha = alpha[1:nspecies])

  #}
  #
  p11 ~ dunif(0.7,1)
  p22 ~ dunif(0.7,1)
  omega[1,1] <- p11
  omega[1,2] <- 1- p11
  omega[2,1] <- 1- p22
  omega[2,2] <- p22
  #omega[1,2] <- 1- omega[1,1]
  #omega[2,1] <- 1- omega[2,2]
  #
  # #omega[1:nspecies, 1:nspecies] <- nimble_omega(p11, p22)
  #
  # r[1:nsite,1:20] <- nimble_INLA(omega[1:nspecies,1:nspecies]) #change 38 to a constant to be specified
  #
  # for(i in 1:nsite){
  #   # for(j in 1:nspecies){
  #   log(lambda_obs[i,1]) <- r[i,5] + r[i,6]*true_cov[i] + r[i,11]+
  #     r[i, 7] + r[i,8]*bias_cov[i]+ r[i,12] - log(1+exp(r[i, 7] + r[i,8]*bias_cov[i]+ r[i,12])) +
  #     r[i, 9] + r[i,10]*det_cov[i] - log(1+exp(r[i, 9] + r[i,10]*det_cov[i]))
  #   # Second species
  #   log(lambda_obs[i,2]) <- r[i,13] + r[i,14]*true_cov[i] + r[i,19]+
  #     r[i, 15] + r[i,16]*bias_cov[i]+ r[i,20] - log(1+exp(r[i, 15] + r[i,16]*bias_cov[i]+ r[i,20])) +
  #     r[i, 17] + r[i,18]*det_cov[i] - log(1+exp(r[i, 17] + r[i,18]*det_cov[i]))
  #
  # }
  #
  # lambda[1:nsite, 1:nspecies] <- lambda_obs[1:nsite, 1:nspecies]
  #
  # #Proportion of lambdas
  #
  #
  # # Proportion for the multinomial distribution
  # for(site.tag in 1:nsite){
  #   for(spe.tag in 1:nspecies){
  #     prop[site.tag,spe.tag] <- (lambda[site.tag, spe.tag])/sum(lambda[site.tag, 1:nspecies])
  #   }
  # }


  # True data
  #for(site.tag in 1:nsite){
    Y[1:nsite] ~ dSpatial(p11,p22, C_obs[1:nsite])
  #}

    prop[1:nsite, 1:nspecies] <- nimblePropTrueClass(p11, p22, Y[1:nsite])

  # Reported species
  for(site.tag in 1:nsite){
    C[site.tag] ~ dcat(prop[site.tag, 1:nspecies])
  }

})


#Data
inla_data <- list(Y=tmp$true_species,
                  C = tmp$error,
                  C_obs = tmp$error)

#Constants
const <- list(nspecies = 2,
              nsite = length(tmp$true_species),
              alpha = c(1,1)
)

# Initial values
idm_inits <- function(){list(p11 = 0.8,
                             p22 = 0.9
  #omega = matrix(c(0.5, 0.5,
                               #               0.5, 0.5), 2,2)
)
}

initsList <- idm_inits()

#Putting all together for the creating compilation
modelInfo <- list(
  code = code,
  constants = const,
  data = inla_data,
  inits = initsList
)

propTrueClass(0.8, 1,tmp$true_species )
#Create the model in nimble
mwtc <- nimbleModel(code,
                    data = inla_data,
                    constants = const,
                    inits = initsList)
#)
#library(igraph)
#plot(mwtc$modelDef$graph)

# Create the model in C
Cmwtc <- compileNimble(mwtc,
                       showCompilerOutput = FALSE) #Have issues compiling


mcmcconf <- configureMCMC(Cmwtc,
                          print=TRUE,
                          useConjugacy=FALSE,
                          monitors = c("omega", "p11", "p22"))

#mcmcconf$removeSamplers(c("omega[1,1:2]","omega[2,1:2]"))
#mcmcconf$addSampler(c("omega[1,1:2]"), "myRW_dirichlet")
# mcmcconf$addSampler(c("omega[2,1:2]"), "myRW_dirichlet")
mcmcconf$removeSamplers(c("p11", "p22"))
mcmcconf$addSampler(c("p11", "p22"), "RW_block",
                    print = TRUE,
                    control = list(scale = 0.5))



Rmcmc <- buildMCMC(mcmcconf)
#enableWAIC =FALSE)

# Compile
cmcmc <- compileNimble(Rmcmc,
                       project = Cmwtc)

# Run the MCMC
#library(pbapply)

mcmc.out <- runMCMC(cmcmc,
                    niter = 100,
                    nburnin = 10 ,
                    inits = initsList,
                    thin = 1,
                    nchains = 1,
                    # thin =100,
                    #setSeed = x,
                    samples=TRUE,
                    samplesAsCodaMCMC = TRUE,
                    summary = TRUE,
                    WAIC = FALSE)
save(mcmc.out, file = paste0("VSEDetectMisclassProbs/VSEDetectMisclass",iter,".RData"))
}



