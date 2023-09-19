#load data
library(INLA)
library(inlabru)
library(dplyr)
library(sf)
library(terra)
bru_safe_sp(force = TRUE)
bru_options_set(control.compute = list(dic = TRUE,
                                       cpo = TRUE,
                                       waic = TRUE))

load("allDataForModelNew.RData")

#Extract parameters needed
covariates = allDataForModels$covariates
mesh = allDataForModels$mesh
boundary = allDataForModels$boundary
dataBeforeClassification = allDataForModels$dataBeforeClassification
spPointsGbif = allDataForModels$spPointsGbif
detdata = allDataForModels$detdata
nspecies = allDataForModels$nspecies
spdeslist = allDataForModels$spdeslist
countDf = allDataForModels$countDf
nbic <- allDataForModels$nbicDf

#Classification probability
omega <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

elev <- covariates[[3]]
elev$alt <- (elev$alt - mean(elev$alt, na.rm = TRUE))/sd(elev$alt, na.rm = TRUE)
f.elev <- function(x, y) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_CRS(elev))
  proj4string(spp) <- fm_CRS(elev)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, elev)
  if (any(is.na(v$elevation))) {
    v$elevation <- inlabru:::bru_fill_missing(elev, spp, v$elevation)
  }
  return(v$elevation)
}

#dist <- raster(distanceTrondheim)
dist <- covariates[[2]]
dist$distanceRaster <- (dist$distanceRaster - mean(dist$distanceRaster, na.rm = TRUE))/sd(dist$distanceRaster, na.rm = TRUE)
f.dist <- function(x, y) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_CRS(dist))
  proj4string(spp) <- fm_CRS(dist)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, dist)
  if (any(is.na(v$layer))) {
    v$layer <- inlabru:::bru_fill_missing(dist, spp, v$layer)
  }
  return(v$layer)
}

prec <- covariates[[1]]
prec$prec01 <- (prec$prec01 - mean(prec$prec01, na.rm = TRUE))/sd(prec$prec01, na.rm = TRUE)
f.prec <- function(x, y) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = CRS(prec))
  proj4string(spp) <- fm_CRS(prec)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, prec)
  if (any(is.na(v$prec01))) {
    v$prec01 <- inlabru:::bru_fill_missing(prec, spp, v$prec01)
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
                resolution=c(3,3), vals=0)
#ras_dom
timeSpent <- rasterize(timeSpentAsSf, ras_dom , "value", update = TRUE) 
  #rasterize(timeSpentAsSf)#, r, "value", background = 0)

distwb <- covariates[[4]]
distwb$distanceRasterWaterBody <- (distwb$distanceRasterWaterBody - mean(distwb$distanceRasterWaterBody, na.rm = TRUE))/sd(distwb$distanceRasterWaterBody, na.rm = TRUE)
f.distwb <- function(x, y) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = CRS(prec))
  proj4string(spp) <- fm_CRS(prec)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, prec)
  if (any(is.na(v$prec01))) {
    v$prec01 <- inlabru:::bru_fill_missing(prec, spp, v$prec01)
  }
  return(v$prec01)
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

cmp1 <- list()
for(i in 1:nspecies){
  cmp1[[i]] <- (paste0("+ betaPO0",i,
                       "(1)+ beta0thin(1)+cov2(main = dist, model='linear')",
                       "+w1",i,"(main = coordinates, model =","spdes[[",i, "]])+",
                       "w2(main = coordinates, model = spde2)+",
                       "cov3", i, "(main = prec, model = 'linear')",
                       "+cov4", i, "(main = distwb, model = 'linear')",
                       "+cov1",i, "(main=elev, model='linear') +beta0det",
                       i,"(1) +timeSpent",i, "(main=timeSpent, model='linear')",
                       "+ betaNina0",i,"(1)+ betaNibc", i, "(1)+betaPA", i, "(1)")
  )
}
cmp <- as.formula(paste0("~ -1",do.call("paste0",cmp1)))

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

spdes <- spdeslist$spdes
spde2 <- spdeslist$spde2

aa <- boundary
#ips <- ipoints(boundary, mesh)
####################################
# Fitting the model
####################################

#}

spdes <- spdeslist$spdes
spde2 <- spdeslist$spde2

aa <- boundary
#ips <- ipoints(boundary, mesh)
#Likelihoods of INLABru
lik1 <- lik3 <- lik4 <- list()
for(i in 1:nspecies){
  lik1[[i]] <- inlabru::like("cp",
                             formula = as.formula(paste0("coordinates ~ beta0",i,"  + cov1",
                                                         i,"  + cov4",
                                                         i," + w1",i,"+cov3",i, "+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)+ beta0det",
                                                         i,"+timeSpent",i,"+fun1(beta0det",i,", timeSpent",i,") + log(weight)")),
                             data = dataBeforeClassification[[i]],
                             #components = cmp,
                             domain = list(coordinates = mesh),
                             #ips = ips
                             samplers = aa
  )
  # lik2[[i]] <- inlabru::like("cp",
  #                            formula = coordinates ~ beta0thin + cov2 + w2,
  #                            data = spPointsGbif,
  #                            #components = cmp,
  #                            domain = list(coordinates = mesh),
  #                            #ips = ips#,
  #                            samplers = aa
  # )
  lik3[[i]] <- inlabru::like("binomial",
                             formula = as.formula(paste0("detdata",i," ~ beta0det",i," + timeSpent",i)),
                             data = detdata[[i]],
                             Ntrials = detdata[[i]]$Ntrials,
                             control.family= list(link = 'logit'),
                             #ips = ips,
                             #components = cmp,
                             domain = list(coordinates = mesh),
                             samplers = aa
  )
  
  lik4[[i]] <- inlabru::like("binomial",
                             formula = as.formula(paste0("detdata",i," ~ beta0",i,"  + cov1",
                                                         i,"+cov3",i, "+cov4",i)),
                             data = detdata1[[i]],
                             #Ntrials = detdata[[i]]$Ntrials,
                             control.family= list(link = 'cloglog'),
                             #ips = ips,
                             #components = cmp,
                             domain = list(coordinates = mesh),
                             samplers = aa
  )
}

elev$cov11 <- elev$cov12  <- elev$alt
dist$cov2 <- dist$distanceRaster
prec$cov31 <- prec$cov32 <- prec$prec01
timeSpent$timeSpent1 <- timeSpent$timeSpent2 <- timeSpent$last
distwb$cov41 <- distwb$cov42 <- distwb$distanceRasterWaterBody

fit5 <- inlabru::bru(cmp,
                     lik1[[1]], lik1[[2]],
                     #unquote(paste0("lik1[[", 1:nspecies, "]]", collapse = ",")),
                     #lik2[[1]],
                     #unquote(paste0("lik3[[", 1:nspecies, "]]", collapse = ",")),
                     lik3[[1]],lik3[[2]],
                     lik4[[1]],lik4[[2]],
                     options = list(control.inla = list(int.strategy = "eb",
                                                        control.vb=list(enable=FALSE),
                                                        cmin = 0.01),
                                    bru_method = list(
                                      #taylor = "pandemic",
                                      # search = "all",
                                      # factor = (1 + sqrt(5)) / 2,
                                      rel_tol = 0.01#,
                                      # max_step = 2,
                                      # lin_opt_method = "onestep"
                                    ), #change to 0.01
                                    bru_max_iter = 1#,
                                   # bru_verbose = 2
                                    )
)

save(fit5, file = "VSEDetectMisclassFit3.RData")

# e.int <- predict(fit5, pixels(mesh, mask = boundary), ~ exp(beta01 + cov31 + cov11+ w11 + beta0thin + cov2 + beta0det1 + timeSpent1))
# 
# ggplot() +
#   gg(e.int) +
#   gg(boundary, alpha = 0) +
#   #gg(nests, shape = "+") +
#   coord_equal()

