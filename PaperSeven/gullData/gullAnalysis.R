#load data
library(INLA)
library(inlabru)
library(dplyr)
library(sf)
library(terra)
library(raster)
library(dplyr)
library(fmesher)
bru_options_set(inla.mode = "experimental")
#bru_safe_sp(force = TRUE)
bru_options_set(control.compute = list(dic = TRUE,
                                       cpo = TRUE,
                                       waic = TRUE))

load("allDataForModelNew.RData")

#Extract parameters needed
BNGproj <- CRS("+proj=robin +datum=WGS84 +units=km")
covariates = allDataForModels$covariates
mesh = allDataForModels$mesh
boundary = allDataForModels$boundary
dataBeforeClassification = lapply(allDataForModels$dataBeforeClassification, function(x){x%>%
    st_as_sf()})
spPointsGbif = allDataForModels$spPointsGbif
detdata = allDataForModels$detdata
nspecies = allDataForModels$nspecies
spdeslist = allDataForModels$spdeslist
countDf = allDataForModels$countDf
nbic <- allDataForModels$nbicDf
occurenceDf <- allDataForModels$occurenceDF

#creating mesh
allData <- do.call("rbind", dataBeforeClassification)
max.edge = diff(range(st_coordinates(allData)[,1])/(3*5))
bound.outer <- diff(range(st_coordinates(allData)[,1])/(5))

#proj4string = CRS(projection(norwaySP))
mesh <- fm_mesh_2d(boundary = boundary,
                   loc = st_coordinates(allData),
                   max.edge  = c(0.5,10)*max.edge,
                   offset = c(max.edge+5,bound.outer),
                   cutoff = 8,
                   #min.angle = 60,
                   crs = BNGproj
)


#df <- pixels(mesh, mask = norwaySP )
#identicalCRS(mesh,boundary)
plot(mesh)
points(dataBeforeClassification[[1]])













#Classification probability
omega <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

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
  v$alt[is.na(v$alt)] <- min(v, na.rm = TRUE)#0#bru_fill_missing(elev, spp, v$alt)
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
  v$distanceRaster[is.na(v$distanceRaster)] <- max(v, na.rm = TRUE)#inlabru:::bru_fill_missing(dist, spp, v$layer)
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
    v$prec01[is.na(v$prec01)] <-max(v, na.rm = TRUE)#inlabru:::bru_fill_missing(prec, spp, v$prec01)
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
timeSpent.spix$layer <- log(timeSpent.spix$layer)# + 0.0001)# - mean(timeSpent.spix$layer))/sd(timeSpent.spix$layer)
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
for(i in 1:2){
  detdata[[i]] <- detdata[[i]]%>%
    #ungroup()%>%
    st_as_sf()
  
  occurenceDf[[i]] <- occurenceDf[[i]]%>%
    #ungroup()%>%
    st_as_sf()
  
  dataBeforeClassification[[i]] <- dataBeforeClassification[[i]][, c("geometry", "weight")]
}
#add covariates
countDf = allDataForModels$countDf
nbic <- allDataForModels$nbicDf

aa <- boundary%>%
  st_as_sf()

cmp1 <- list()

for(i in 1:nspecies){
  
  cmp1[[i]] <- paste0("+betaNina0",i,"(1) + cov1",i, "(f.elev(st_coordinates(.data.)), model = 'linear')+cov3", i, "(f.prec(st_coordinates(.data.)), model = 'linear')+cov4", i,"(f.distwb(st_coordinates(.data.)), model = 'linear')+betaNibc0",i, "(1)+ latEffect",i, "(main = decimalLatitude, model = 'iid')+yearEffect", i,"(main = year, model = 'iid')","+beta0det",i,"(1) + timeSpent",i,"(f.time(st_coordinates(.data.)), model = 'linear')+ betaPo0",i,"(1) + w1", i, "(main = geometry, model =","spdes[[",i, "]]) + betaPA0", i, "(1) + cov2",i,"(f.dist(st_coordinates(.data.)), model = 'linear') + w2(main = geometry, model =spde2)+ weight(f.weights(.data.,st_coordinates(.data.)), model = 'offset')")
}
cmp <- as.formula(paste0("~ -1",do.call("paste0",cmp1)))

#cmp <- as.formula(~-1 + betaNina01(1) + betaNina02(1) )
####################################
# Fitting the model
####################################
fun <- function(x,y,z){
  -log(1+exp((x+y+z)))
}

fun1 <- function(x,y){
  -log(1+exp((x+y)))
}

#if(nspecies == 2){
functionOmega <-  function(a, d, b, e, c, f,g, h, x){
  ret <- log(omega[1,x]*plogis(a+b+c+g) + omega[2,x]*plogis(d+e+f+h))
  return(ret)
}
#}

range <- diff(range(mesh$loc[,2])) 
spdeEological <- list()
#for(i in 1: nspecies){
spdeEological[[1]] <- inla.spde2.pcmatern(mesh = mesh,
                                          # PC-prior on range: P(practic.range < 0.05) = 0.01
                                          prior.range = c(10, 0.1),
                                          # PC-prior on sigma: P(sigma > 0.8) = 0.01
                                          prior.sigma = c(1, 0.1))

spdeEological[[2]] <- inla.spde2.pcmatern(mesh = mesh,
                                          # PC-prior on range: P(practic.range > range) = 0.95
                                          prior.range = c(10, 0.1),
                                          # PC-prior on sigma: P(sigma > 3) = 0.01
                                          prior.sigma = c(1, 0.1))
#}

#SPDEs for the thinning
spdeSampling <- inla.spde2.pcmatern(mesh = mesh,
                                    # PC-prior on range: P(practic.range < 0max.edge+3) = 0.01
                                    prior.range = c(10, 0.1),
                                    # PC-prior on sigma: P(sigma > 0.8) = 0.01
                                    prior.sigma = c(1, 0.1))
spdeslist <- list(spdes = spdeEological, spde2 = spdeSampling)

spdes <- spdeslist$spdes
spde2 <- spdeslist$spde2

aa <- boundary%>%
  st_as_sf()


