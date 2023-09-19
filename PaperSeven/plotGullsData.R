# Load packages
load("gulldata/allDataForModelNew.RData")
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

#source("DataGeneration.R")
library(sf)
library(terra)
library(raster)
library(inlabru)
library(ggplot2)
library(INLA)
library(dplyr)
library(fmesher)
#Extract parameters needed
covariates = allDataForModels$covariates
mesh = allDataForModels$mesh
boundary = allDataForModels$boundary
detdata = allDataForModels$detdata
nspecies = allDataForModels$nspecies
spdeslist = allDataForModels$spdeslist
countDf = allDataForModels$countDf
nbic <- allDataForModels$nbicDf
occurenceDf <- allDataForModels$occurenceDF
nspecies <- 2

#mesh used to fit data
BNGproj <- CRS("+proj=robin +datum=WGS84 +units=km")
dataBeforeClassification = lapply(allDataForModels$dataBeforeClassification, function(x){x%>%
    st_as_sf()})
#creating mesh
allData <- do.call("rbind", dataBeforeClassification)
max.edge = diff(range(st_coordinates(allData)[,1])/(3*5))
bound.outer <- diff(range(st_coordinates(allData)[,1])/(5))

#proj4string = CRS(projection(norwaySP))
mesh <- fm_mesh_2d(boundary = boundary,
                   loc = st_coordinates(allData),
                   max.edge  = c(20,30)*max.edge,
                   offset = c(max.edge+5,bound.outer),
                   cutoff = 8,
                   #min.angle = 60,
                   crs = BNGproj
)

elev <- covariates[[3]]
elev$alt <- elev$alt/1000#- mean(elev$alt, na.rm = TRUE))/sd(elev$alt, na.rm = TRUE)
#elev$alt <- (elev$alt - mean(elev$alt, na.rm = TRUE))/sd(elev$alt, na.rm = TRUE)
f.elev <- function(xy) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  xy1 <- xy#st_coordinates(xy)
  spp <- SpatialPoints(data.frame(x = xy1[,1], y = xy1[, 2]),
                       proj4string = fm_CRS(elev)
  )
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, elev)
  #if (any(is.na(v$alt))) {
  v$alt[is.na(v$alt)] <- 0#0#bru_fill_missing(elev, spp, v$alt)
  #}
  return(v$alt)
}


#dist <- raster(distanceTrondheim)
dist <- covariates[[2]]
dist$distanceRaster <- dist$distanceRaster/1000
#dist$distanceRaster <- (dist$distanceRaster- mean(dist$distanceRaster, na.rm = TRUE))/sd(dist$distanceRaster, na.rm = TRUE)
f.dist <- function(xy) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  xy1 <- xy#st_coordinates(xy)
  spp <- SpatialPoints(data.frame(x = xy1[,1], y = xy1[,2]), proj4string = fm_CRS(dist))
  #proj4string(spp) <- fm_CRS(dist)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, dist)
  #if (any(is.na(v$distanceRaster))) {
  v$distanceRaster[is.na(v$distanceRaster)] <- 0#inlabru:::bru_fill_missing(dist, spp, v$layer)
  #}
  return(v$distanceRaster)
}

prec <- covariates[[1]]
prec$prec01 <- prec$prec01/1000# - mean(prec$prec01, na.rm = TRUE))/sd(prec$prec01, na.rm = TRUE)
#prec$prec01 <- (prec$prec01 - mean(prec$prec01, na.rm = TRUE))/sd(prec$prec01, na.rm = TRUE)

f.prec <- function(xy) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  xy1 <- xy#st_coordinates(xy)
  spp <- SpatialPoints(data.frame(x = xy1[,1], y = xy1[,2]), proj4string = fm_CRS(prec))
  #proj4string(spp) <- fm_CRS(prec)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, prec)
  if (any(is.na(v$prec01))) {
    v$prec01[is.na(v$prec01)] <- 0.16#inlabru:::bru_fill_missing(prec, spp, v$prec01)
  }
  return(v$prec01)
}

timeSpent1 <- detdata[[1]]%>%
  st_as_sf()
coords <- st_coordinates(timeSpent1)
timeSpentNew <- data.frame(x = coords[,1],
                           y = coords[,2],
                           value = timeSpent1$timeSpent1)

timeSpentAsSf <- timeSpentNew%>%
  st_as_sf(., coords = c("x", "y"), crs = CRS("+proj=robin +datum=WGS84 +units=km"))

ext <- c(min(mesh$loc[,1]), max(mesh$loc[,1]), min(mesh$loc[,2]), max(mesh$loc[,2]))

ras_dom <-raster(xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4],
                 crs="+proj=robin +datum=WGS84 +units=km",
                 resolution=c(3,3), vals=1)
#ras_dom
timeSpent <- rasterize(timeSpentAsSf, ras_dom , "value", update = TRUE) 
timeSpent.spix <- as(timeSpent,"SpatialPixelsDataFrame")
timeSpent.spix$layer <- log(timeSpent.spix$layer)#/1000# + 0.0001)# - mean(timeSpent.spix$layer))/sd(timeSpent.spix$layer)
#timeSpent.spix$layer <- (timeSpent.spix$layer - mean(timeSpent.spix$layer))/sd(timeSpent.spix$layer)
f.time <- function(xy){
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  xy1 <- xy#st_coordinates(xy)
  spp <- SpatialPoints(data.frame(x = xy1[,1], y = xy1[,2]), proj4string = fm_CRS(dist))
  #proj4string(spp) <- fm_CRS(prec)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, timeSpent.spix)
  #if (any(is.na(v$layer))) {
  v$layer[is.na(v$layer)] <- 0#inlabru:::bru_fill_missing(prec, spp, v$prec01)
  #}
  return(v$layer)
}

#rasterize(timeSpentAsSf)#, r, "value", background = 0)

distwb <- covariates[[4]]
distwb$distanceRasterWaterBody <- distwb$distanceRasterWaterBody/1000 #- mean(distwb$distanceRasterWaterBody, na.rm = TRUE))/sd(distwb$distanceRasterWaterBody, na.rm = TRUE)
#distwb$distanceRasterWaterBody <- (distwb$distanceRasterWaterBody - mean(distwb$distanceRasterWaterBody))#/sd(distwb$distanceRasterWaterBody)
f.distwb <- function(xy) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  xy1 <- xy#st_coordinates(xy)
  spp <- SpatialPoints(data.frame(x = xy1[,1], y = xy1[,2]), proj4string = fm_CRS(distwb))
  #spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = CRS(prec))
  #proj4string(spp) <- fm_CRS(prec)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, distwb)
  #if (any(is.na(v$distanceRasterWaterBody))) {
  v$distanceRasterWaterBody[is.na(v$distanceRasterWaterBody)] <- 0#inlabru:::bru_fill_missing(prec, spp, v$prec01)
  #}
  
  
  return(v$distanceRasterWaterBody)
}

f.weights <- function(dd, xy) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  xy1 <- xy #st_coordinates(xy)
  ddDf <- as.data.frame(dd)%>%
    dplyr::mutate(X = st_coordinates(geometry)[,1],
                  Y = st_coordinates(geometry)[,2])
  weight <- log(ddDf[rowSums(ddDf[,3:4] == xy1[,1:2])==2, 2])
  
  #if (any(is.na(weight))) {
  weight[is.na(weight)] <- 0#inlabru:::bru_fill_missing(prec, spp, v$prec01)
  #}
  
  
  return(weight)
}
# f.timeSpent <- function(x, y) {
#   # turn coordinates into SpatialPoints object:
#   # with the appropriate coordinate reference system (CRS)
#   spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(timeSpent))
#   proj4string(spp) <- fm_sp_get_crs(timeSpent)
#   # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
#   v <- over(spp, dist)
#   if (any(is.na(v$layer))) {
#     v$layer <- inlabru:::bru_fill_missing(dist, spp, v$layer)
#   }
#   return(v$layer)
# }
# for(i in 1:2){
#   detdata[[i]] <- detdata[[i]]%>%
#     #ungroup()%>%
#     st_as_sf()
#   
#   occurenceDf[[i]] <- occurenceDf[[i]]%>%
#     #ungroup()%>%
#     st_as_sf()
#   
#   dataBeforeClassification[[i]] <- dataBeforeClassification[[i]][, c("geometry", "weight")]
# }
#add covariates
countDf = allDataForModels$countDf
nbic <- allDataForModels$nbicDf

##############################
# Naive model Fit
#############################

# load("gullData/naiveModelFit3.RData")
# tmp <- fit1$marginals.fixed
# betaPo0 <- betaPA0 <- betaNina0 <- betaNibc0 <- beta1 <- beta2 <- vector("numeric", 2)
# for(i in 1:nspecies){
#   betaPo0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPo0",i)][[1]])
#   betaPA0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPA0",i)][[1]])
#   betaNina0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNina0",i)][[1]])
#   betaNibc0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNibc0",i)][[1]])
#   beta1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]])
#   beta2[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
# }
# 
# gamma0 <- gamma1 <- c(NA, NA)
# alpha0 <- alpha1 <- NA
# 
# ## The locations where information is needed
# tmpHyper <- fit1$marginals.hyperpar
# rangeW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w11`)
# stdW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w11`)
# rangeW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w12`)
# stdW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w12`)
# rangeW2 = NA
# stdW2 = NA
# mld <- fit1$mlik[1,1]
# dic <- fit1$dic$dic
# waic <- fit1$waic$waic
# fittedNaive <- data.frame(mld = mld,
#                          dic = dic,
#                          waic = waic,
#                          alpha0 = alpha0,
#                          alpha1 = alpha1,
#                          betaPo01 = betaPo0[1],
#                          betaPo02 = betaPo0[2],
#                          betaPA01 = betaPA0[1],
#                          betaPA02 = betaPA0[2],
#                          betaNina01 = betaNina0[1],
#                          betaNina02 = betaNina0[2],
#                          betaNibc01 = betaNibc0[1],
#                          betaNibc02 = betaNibc0[2],
#                          beta11 = beta1[1],
#                          beta12 = beta1[2],
#                          beta21 = beta2[1],
#                          beta22 = beta2[2],
#                          gamma01 = gamma0[1],
#                          gamma02 = gamma0[2],
#                          gamma11 = gamma1[1],
#                          gamma12 = gamma1[2],
#                          rangeW11 = rangeW11,
#                          stdW11 = stdW11,
#                          rangeW12 = rangeW12,
#                          stdW12 = stdW12,
#                          rangeW2 = rangeW2,
#                          stdW2 = stdW2)
# 
# 
# ########
# # VSE
# ##########
# load("gullData/samplingEffortFit3.RData")
# gamma0 <- gamma1 <- c(NA, NA)
# tmp <- fit2$marginals.fixed
# alpha0 <- NA#c(INLA::inla.emarginal(function(x) x,tmp[paste0("beta0thin")][[1]]))
# alpha1 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("cov2")][[1]]))
# betaPo0 <- betaPA0 <- betaNina0 <- betaNibc0 <- beta1 <- beta2 <- vector("numeric", 2)
# for(i in 1:nspecies){
#   betaPo0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPo0",i)][[1]])
#   betaPA0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPA0",i)][[1]])
#   betaNina0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNina0",i)][[1]])
#   betaNibc0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNibc0",i)][[1]])
#   beta1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]])
#   beta2[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
# }
# gamma0 <- gamma <- c(NA, NA)
# 
# 
# ## The locations where information is needed
# tmpHyper <- fit2$marginals.hyperpar
# rangeW2 <- NA#INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w2`)
# stdW2 <- NA #INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w2`)
# rangeW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w11`)
# stdW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w11`)
# rangeW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w12`)
# stdW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w12`)
# mld <- fit2$mlik[1,1]
# dic <- fit2$dic$dic
# waic <- fit2$waic$waic
# fittedVSE <- data.frame(mld = mld,
#                         dic = dic,
#                         waic = waic,
#                         alpha0 = alpha0,
#                         alpha1 = alpha1,
#                         betaPo01 = betaPo0[1],
#                         betaPo02 = betaPo0[2],
#                         betaPA01 = betaPA0[1],
#                         betaPA02 = betaPA0[2],
#                         betaNina01 = betaNina0[1],
#                         betaNina02 = betaNina0[2],
#                         betaNibc01 = betaNibc0[1],
#                         betaNibc02 = betaNibc0[2],
#                         beta11 = beta1[1],
#                         beta12 = beta1[2],
#                         beta21 = beta2[1],
#                         beta22 = beta2[2],
#                         gamma01 = gamma0[1],
#                         gamma02 = gamma0[2],
#                         gamma11 = gamma1[1],
#                         gamma12 = gamma1[2],
#                         rangeW11 = rangeW11,
#                         stdW11 = stdW11,
#                         rangeW12 = rangeW12,
#                         stdW12 = stdW12,
#                         rangeW2 = rangeW2,
#                         stdW2 = stdW2)
# 
# ##########
# #   Detect Probability
# #############
# load("gullData/detectionProbabilityFit3.RData")
# 
# tmp <- fit3$marginals.fixed
# alpha0 <- alpha1 <- NA
# beta0 <- beta1 <- beta2 <- vector("numeric", 2)
# gamma0 <- gamma1 <- c(NA, NA)
# betaPo0 <- betaPA0 <- betaNina0 <- betaNibc0 <- beta1 <- beta2 <- vector("numeric", 2)
# for(i in 1:nspecies){
#   betaPo0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPo0",i)][[1]])
#   betaPA0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPA0",i)][[1]])
#   betaNina0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNina0",i)][[1]])
#   betaNibc0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNibc0",i)][[1]])
#   beta1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]])
#   beta2[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
#   gamma0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("beta0det",i)][[1]])
#   gamma1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
#   }
# 
# 
# ## The locations where information is needed
# tmpHyper <- fit3$marginals.hyperpar
# rangeW2 <- NA
# stdW2 <- NA
# rangeW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w11`)
# stdW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w11`)
# rangeW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w12`)
# stdW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w12`)
# mld <- fit3$mlik[1,1]
# dic <- fit3$dic$dic
# waic <- fit3$waic$waic
# fittedDetect <- data.frame(mld = mld,
#                            dic = dic,
#                            waic = waic,
#                            alpha0 = alpha0,
#                            alpha1 = alpha1,
#                            betaPo01 = betaPo0[1],
#                            betaPo02 = betaPo0[2],
#                            betaPA01 = betaPA0[1],
#                            betaPA02 = betaPA0[2],
#                            betaNina01 = betaNina0[1],
#                            betaNina02 = betaNina0[2],
#                            betaNibc01 = betaNibc0[1],
#                            betaNibc02 = betaNibc0[2],
#                            beta11 = beta1[1],
#                            beta12 = beta1[2],
#                            beta21 = beta2[1],
#                            beta22 = beta2[2],
#                            gamma01 = gamma0[1],
#                            gamma02 = gamma0[2],
#                            gamma11 = gamma1[1],
#                            gamma12 = gamma1[2],
#                            rangeW11 = rangeW11,
#                            stdW11 = stdW11,
#                            rangeW12 = rangeW12,
#                            stdW12 = stdW12,
#                            rangeW2 = rangeW2,
#                            stdW2 = stdW2)
# ################
# # VSEDetect
# ###############
# load("gullData/VSEDetectFit3.RData")
# tmp <- fit4$marginals.fixed
# alpha0 <- NA#c(INLA::inla.emarginal(function(x) x,tmp[paste0("beta0thin")][[1]]))
# alpha1 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("cov2")][[1]]))
# betaPo0 <- betaPA0 <- betaNina0 <- betaNibc0 <- beta1 <- beta2 <- vector("numeric", 2)
# for(i in 1:nspecies){
#   betaPo0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPo0",i)][[1]])
#   betaPA0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPA0",i)][[1]])
#   betaNina0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNina0",i)][[1]])
#   betaNibc0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNibc0",i)][[1]])
#   beta1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]])
#   beta2[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
#   gamma0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("beta0det",i)][[1]])
#   gamma1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
# }
# 
# 
# ## The locations where information is needed
# tmpHyper <- fit3$marginals.hyperpar
# rangeW2 <- NA
# stdW2 <- NA
# rangeW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w11`)
# stdW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w11`)
# rangeW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w12`)
# stdW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w12`)
# mld <- fit3$mlik[1,1]
# dic <- fit3$dic$dic
# waic <- fit3$waic$waic
# fittedVSEDetect <- data.frame(mld = mld,
#                               dic = dic,
#                               waic = waic,
#                               alpha0 = alpha0,
#                               alpha1 = alpha1,
#                               betaPo01 = betaPo0[1],
#                               betaPo02 = betaPo0[2],
#                               betaPA01 = betaPA0[1],
#                               betaPA02 = betaPA0[2],
#                               betaNina01 = betaNina0[1],
#                               betaNina02 = betaNina0[2],
#                               betaNibc01 = betaNibc0[1],
#                               betaNibc02 = betaNibc0[2],
#                               beta11 = beta1[1],
#                               beta12 = beta1[2],
#                               beta21 = beta2[1],
#                               beta22 = beta2[2],
#                               gamma01 = gamma0[1],
#                               gamma02 = gamma0[2],
#                               gamma11 = gamma1[1],
#                               gamma12 = gamma1[2],
#                               rangeW11 = rangeW11,
#                               stdW11 = stdW11,
#                               rangeW12 = rangeW12,
#                               stdW12 = stdW12,
#                               rangeW2 = rangeW2,
#                               stdW2 = stdW2)
# 
# ################
# # VSEDetectMisclass
# ###############
# load("gullData/VSEDetectMisclassFit3.RData")
# tmp <- fit5$marginals.fixed
# alpha0 <-NA #c(INLA::inla.emarginal(function(x) x,tmp[paste0("beta0thin")][[1]]))
# alpha1 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("cov2")][[1]]))
# betaPo0 <- betaPA0 <- betaNina0 <- betaNibc0 <- beta1 <- beta2 <- vector("numeric", 2)
# for(i in 1:nspecies){
#   betaPo0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPo0",i)][[1]])
#   betaPA0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaPA0",i)][[1]])
#   betaNina0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNina0",i)][[1]])
#   betaNibc0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("betaNibc0",i)][[1]])
#   beta1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]])
#   beta2[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
#   gamma0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("beta0det",i)][[1]])
#   gamma1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
# }
# 
# 
# ## The locations where information is needed
# tmpHyper <- fit5$marginals.hyperpar
# rangeW2 <- NA#INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w2`)
# stdW2 <- NA#INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w2`)
# rangeW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w11`)
# stdW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w11`)
# rangeW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w12`)
# stdW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w12`)
# mld <- fit5$mlik[1,1]
# dic <- fit5$dic$dic
# waic <- fit5$waic$waic
# fittedVSEDetectMisclass <- data.frame(mld = mld,
#                                       dic = dic,
#                                       waic = waic,
#                                       alpha0 = alpha0,
#                                       alpha1 = alpha1,
#                                       betaPo01 = betaPo0[1],
#                                       betaPo02 = betaPo0[2],
#                                       betaPA01 = betaPA0[1],
#                                       betaPA02 = betaPA0[2],
#                                       betaNina01 = betaNina0[1],
#                                       betaNina02 = betaNina0[2],
#                                       betaNibc01 = betaNibc0[1],
#                                       betaNibc02 = betaNibc0[2],
#                                       beta11 = beta1[1],
#                                       beta12 = beta1[2],
#                                       beta21 = beta2[1],
#                                       beta22 = beta2[2],
#                                       gamma01 = gamma0[1],
#                                       gamma02 = gamma0[2],
#                                       gamma11 = gamma1[1],
#                                       gamma12 = gamma1[2],
#                                       rangeW11 = rangeW11,
#                                       stdW11 = stdW11,
#                                       rangeW12 = rangeW12,
#                                       stdW12 = stdW12,
#                                       rangeW2 = rangeW2,
#                                       stdW2 = stdW2)
# 
# # Put all data together
# allResults <- rbind(fittedNaive,
#                     fittedVSE,
#                     fittedDetect,
#                     fittedVSEDetect,
#                     fittedVSEDetectMisclass)
#write.csv(allResults, file = "results/gullDataResultsNew")


#Plots for Ecological process pars
load("gullData/naiveModelFit31.RData")
load("gullData/samplingEffortFit31.RData")
load("gullData/detectionProbabilityFit31.RData")
load("gullData/VSEDetectFit31.RData")
#load("gullData/VSEDetectMisclassFit3.RData")

#all fitted models
allModelsFitted <- list(fit1,
                        fit2,
                        fit3,
                        fit4,
                        fit4 #change to fit5
                        )

estimates <- lapply(allModelsFitted, function(x){
  #samples <- INLA::inla.posterior.sample(10000, x)
  #betaNina01 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$betaNina01)), "betaNina01")
  #betaNina02 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$betaNina02)), "betaNina02")
  #samples$latent[grepl(paste0("betaNibc01"), rownames(samples$latent)),]
  betaNibc01 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$betaNibc01)), "betaNibc01")
  betaNibc02 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$betaNibc02)), "betaNibc02")
  betaPo01 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$betaPo01)), "betaPo01")
  betaPo02 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$betaPo02)), "betaPo02")
  betaPA01 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$betaPA01)), "betaPA01")
  betaPA02 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$betaPA02)), "betaNibc02")
  
  elev1 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$cov11)), "elev1")
  elev2 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$cov12)), "elev2")
  prec1 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$cov31)), "prec1")
  prec2 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$cov32)), "prec2")
  distwb1 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$cov41)), "distwb1")
  distwb2 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$cov42)), "distwb2")
  errorIndicator <- inherits(try( dist1 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$cov21)), "dist1"), 
                                  silent = TRUE),
                             "try-error")
  errorIndicator1 <- inherits(try( dist2 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$cov22)), "dist2"), 
                                  silent = TRUE),
                             "try-error")
  if(errorIndicator) dist1 <- c(rep(NA, 7), "dist1")
  if(errorIndicator1) dist2 <- c(rep(NA, 7), "dist2")
  errorIndicator2 <- inherits(try( timeSpent1 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$timeSpent1)), "timeSpent1"), 
                                  silent = TRUE),
                             "try-error")
  errorIndicator3 <- inherits(try( timeSpent2 <- c(unlist(INLA::inla.zmarginal(x$marginals.fixed$timeSpent2)), "timeSpent2"), 
                                   silent = TRUE),
                              "try-error")
  if(errorIndicator2) timeSpent1 <- c(rep(NA, 7), "timeSpent1")
  if(errorIndicator3) timeSpent2 <- c(rep(NA, 7), "timeSpent2")
  # spde component for lambda
  rho11 <- c(unlist(INLA::inla.zmarginal(x$marginals.hyperpar$`Range for w11`)), "rho11")
  rho12 <- c(unlist(INLA::inla.zmarginal(x$marginals.hyperpar$`Range for w11`)), "rho12")
  sd11 <- c(unlist(INLA::inla.zmarginal(x$marginals.hyperpar$`Stdev for w11`)), "sd11")
  sd12 <- c(unlist(INLA::inla.zmarginal(x$marginals.hyperpar$`Stdev for w12`)), "sd12") 

  ret <- rbind( betaNibc01, betaNibc02,
               betaPo01 ,betaPo02 ,betaPA01 ,betaPA01,
               elev1, elev2, prec1, prec2,
               distwb1, distwb2,dist1, dist2, timeSpent1, timeSpent2,
               rho11, rho12, sd11, sd12
               )
  })%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  dplyr::mutate(model = rep(c("naive", "USE", "Detect", "USEDetect", "USEDetMis"), each = 20))

colnames(estimates)[8] <- "parameter"

#save results
write.csv(estimates, file = "results/gullDataResultsNew.csv")

#### Ecological par estimates
 library(readr)
estimates <- read_csv("results/gullDataResultsNew.csv")
fig <- estimates%>%
  dplyr::filter(!model%in%c("USEDetMis"))%>%
  dplyr::filter(parameter %in% c("elev1", "elev2", "prec1", "prec2", "distwb1", "distwb2", "dist1", "dist2", "timeSpent1", "timeSpent2"))%>%
  ggplot(., mapping = aes(x= parameter, y = as.numeric(quant0.5), col = model))+
  #geom_point(position = position_dodge(width = 0.9))+
  geom_pointrange(aes(ymin = as.numeric(quant0.025), ymax = as.numeric(quant0.975)), position = position_dodge(width = 0.9))+
  #ylim(c(-6.1,0))+
  coord_flip()+
  theme_bw()+
  ylab(expression(paste("Posterior Median"%+-%"95% CI")))+
  xlab("Covariate Effect")+
  theme(legend.position = "bottom")

ggsave("results/estimatedParsGulls.png",
       plot = fig,
       width = 10,
       height = 7,
       units = "in")

estimates%>%
  dplyr::filter(parameter %in% c("distwb1"))

#Model assessment checks
modelAssessment <- lapply(allModelsFitted, function(x){
  waic <- x$waic$waic
  cpoVals <- log(x$cpo$cpo)
  cpo <- sum((cpoVals[!is.infinite(cpoVals)]))
  mlike <- x$mlik[1,1]
  ret <- c(waic, cpo, mlike)
  return(ret)
})%>%
  do.call("rbind", .)%>%
  as.data.frame()
colnames(modelAssessment) <- c("WAIC", "lpml", "marg_like")
modelAssessment$WAIC <-  modelAssessment$WAIC - min(modelAssessment$WAIC) 
modelAssessment$lpml <- max(modelAssessment$lpml) - modelAssessment$lpml 

write.csv(modelAssessment, file = "results/gullDataModelAssessment.csv")

aa <- boundary%>%
  st_as_sf()
# Plot predictions
# Get predicted values
predNaive <- predict(
  fit1,
  fm_pixels(mesh, mask = aa, format = "sf"),
  ~ data.frame(
    lambda1 = (betaPo01  + cov11 + w11+ cov31+cov41),
    lambda2 =(betaPo02  + cov12 + w12+cov32+cov42),
    trueLam1 = exp(betaPo01  +  cov31 + cov11 + w11+cov41 + betaPA01+ betaNibc01),
    trueLam2 = exp(betaPo02  + cov12+  cov32 + w12+cov42 + betaPA02 + betaNibc02)
  )
)

gg1 <- ggplot() +
  gg(tidyr::pivot_longer(predNaive$trueLam1,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  ggtitle("Naive - herring gulls")

gg2 <- ggplot() +
  gg(tidyr::pivot_longer(predNaive$trueLam2,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  ggtitle("Naive - common gulls")

#Obtain standard deviation estimates
gg1SD <- ggplot() +
  gg(tidyr::pivot_longer(predNaive$trueLam1,
                         c(sd),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1, limits = c(0,2.5))+
  #facet_wrap(~quantile)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_bw()+
  xlab("")+
  ylab("")+
  ggtitle("Naive - L. argentatus")

gg2SD <- ggplot() +
  gg(tidyr::pivot_longer(predNaive$trueLam2,
                         c(sd),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_bw()+
  xlab("")+
  ylab("")+
  ggtitle("Naive - L. canus")

predVSE <- predict(
  fit2,
  fm_pixels(mesh, mask = aa, format = "sf"),
  ~ data.frame(
    #lambda1 = exp(beta01  + cov11 + cov31 + w11+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
   # lambda2 =exp(beta02  + cov12 +  cov32 + w12+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
    trueLam1 = exp(betaPo01  +  cov31 + cov11 + w11+cov41 + betaPA01+ betaNibc01),
    trueLam2 = exp(betaPo02  + cov12+  cov32 + w12+cov42 + betaPA02 + betaNibc02)
  )
)

#predVSENew <- predVSE
#predVSENew$trueLam1$mean <-  predVSE$trueLam1$median - predNaive$trueLam1$median
#predVSENew$trueLam2$mean <-  predVSE$trueLam2$median - predNaive$trueLam2$median
gg3 <- ggplot() +
  gg(tidyr::pivot_longer(predVSE$trueLam1,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1, limits = c(0,0.5))+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5))+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  ggtitle("C) USE")

gg4 <- ggplot() +
  gg(tidyr::pivot_longer(predVSE$trueLam2,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile")+
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1, limits = c(0,0.5))+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5))+
  ggtitle("D) USE")


predDetect <- predict(
  fit3,
  fm_pixels(mesh, mask = aa, format = "sf"),
  ~ data.frame(
    #lambda1 = exp(beta01  + cov11 + cov31 + w11+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
    # lambda2 =exp(beta02  + cov12 +  cov32 + w12+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
    trueLam1 = exp(betaPo01  +  cov31 + cov11 + w11+cov41 + betaPA01+ betaNibc01),
    trueLam2 = exp(betaPo02  + cov12+  cov32 + w12+cov42 + betaPA02 + betaNibc02)
  )
)

#predVSENew <- predDetect
#predVSENew$trueLam1$mean <-  predDetect$trueLam1$median - predNaive$trueLam1$median
#predVSENew$trueLam2$mean <-  predDetect$trueLam2$median - predNaive$trueLam2$median

gg5 <- ggplot() +
  gg(tidyr::pivot_longer(predDetect$trueLam1,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  ggtitle("Detect - herring gulls")

gg6 <- ggplot() +
  gg(tidyr::pivot_longer(predDetect$trueLam2,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  ggtitle("Detect - common gulls")


predVSEDetect <- predict(
  fit4,
  fm_pixels(mesh, mask = aa, format = "sf"),
  ~ data.frame(
    #lambda1 = exp(beta01  + cov11 + cov31 + w11+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
    # lambda2 =exp(beta02  + cov12 +  cov32 + w12+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
    trueLam1 = exp(betaPo01  +  cov31 + cov11 + w11+cov41 + betaPA01+ betaNibc01),
    trueLam2 = exp(betaPo02  + cov12+  cov32 + w12+cov42 + betaPA02 + betaNibc02)
  )
)


#predVSENew <- predVSEDetect
#predVSENew$trueLam1$mean <-  predVSEDetect$trueLam1$median - predNaive$trueLam1$median
#predVSENew$trueLam2$mean <-  predVSEDetect$trueLam2$median - predNaive$trueLam2$median

gg7 <- ggplot() +
  gg(tidyr::pivot_longer(predVSEDetect$trueLam1,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  ggtitle("USEDetect - herring gulls")

gg8 <- ggplot() +
  gg(tidyr::pivot_longer(predVSEDetect$trueLam2,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile")+
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  ggtitle("USEDetect - common gulls")

#load("gullData/VSEDetectMisclassFit2.RData")
predVSEDetectMisclass <- predict(
  fit5,
  fm_pixels(mesh, mask = aa, format = "sf"),
  ~ data.frame(
    #lambda1 = exp(beta01  + cov11 + cov31 + w11+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
    # lambda2 =exp(beta02  + cov12 +  cov32 + w12+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
    trueLam1 = exp(betaPo01  +  cov31 + cov11 + w11+cov41 + betaPA01+ betaNibc01),
    trueLam2 = exp(betaPo02  + cov12+  cov32 + w12+cov42 + betaPA02 + betaNibc02)
  )
)


#predVSENew <- predVSEDetectMisclass 
#predVSENew$trueLam1$mean <-  predVSEDetectMisclass$trueLam1$median - predNaive$trueLam1$median
#predVSENew$trueLam2$mean <-  predVSEDetectMisclass$trueLam2$median - predNaive$trueLam2$median

gg9 <- ggplot() +
  gg(tidyr::pivot_longer(predVSEDetectMisclass$trueLam1,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5))+
  ggtitle("I) VSEDetectMisclass")

gg10 <- ggplot() +
  gg(tidyr::pivot_longer(predVSEDetectMisclass$trueLam2,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile")+
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5))+
  ggtitle("J) VSEDetectMisclass")


fig <- ggpubr::ggarrange(gg1, gg5,
                  #gg5,gg7,
                  gg7, gg2,
                  gg6, gg8,
                 #gg8, gg10,
                  ncol =3,
                  nrow = 2#,
                 #common.legend = TRUE,
                 #legend = "bottom"
                 )

ggsave("results/predictedIntensity.png",
       plot = fig,
       width = 10,
       height = 7,
       units = "in")


# Predictions from the best model
predVSE <- predict(
  fit2,
  fm_pixels(mesh, mask = aa, format = "sf"),
  ~ data.frame(
    #lambda1 = exp(beta01  + cov11 + cov31 + w11+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
    # lambda2 =exp(beta02  + cov12 +  cov32 + w12+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)),
    trueLam1 = exp(betaPo01  +  cov31 + cov11 + w11+cov41 + betaPA01+ betaNibc01),
    trueLam2 = exp(betaPo02  + cov12+  cov32 + w12+cov42 + betaPA02 + betaNibc02)
  )
)

#predVSENew <- predVSE
#predVSENew$trueLam1$mean <-  predVSE$trueLam1$median - predNaive$trueLam1$median
#predVSENew$trueLam2$mean <-  predVSE$trueLam2$median - predNaive$trueLam2$median
gg3 <- ggplot() +
  gg(tidyr::pivot_longer(predVSE$trueLam1,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile") +
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1, limits = c(0,0.5))+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  ggtitle("Larus argentatus")

gg4 <- ggplot() +
  gg(tidyr::pivot_longer(predVSE$trueLam2,
                         c('q0.5'),
                         names_to = "quantile", values_to = "value"
  ),
  aes(fill = value),
  geom = "tile")+
  gg(aa, alpha = 0) + #boundary
  scale_fill_viridis_c(direction = -1, limits = c(0,0.5))+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_sf()+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.1, 'cm'),
        legend.text = element_text(size =5),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  ggtitle("Larus canus")

#save plot
fig <- ggpubr::ggarrange(gg3, gg4,
                         ncol =2,
                         nrow = 1#,
                         #common.legend = TRUE#,
                         #legend = "bottom"
)

ggsave("results/predictedIntensityBestModel.png",
       plot = fig,
       width = 10,
       height = 5,
       units = "in")


