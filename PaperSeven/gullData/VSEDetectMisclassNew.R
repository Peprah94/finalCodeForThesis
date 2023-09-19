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
                   max.edge  = c(20,30)*max.edge,
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
                 resolution=c(3,3), vals=0)
#ras_dom
timeSpent <- rasterize(timeSpentAsSf, ras_dom , "value", update = TRUE) 
timeSpent.spix <- as(timeSpent,"SpatialPixelsDataFrame")
timeSpent.spix$layer <- timeSpent.spix$layer/1000# + 0.0001)# - mean(timeSpent.spix$layer))/sd(timeSpent.spix$layer)
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

# f.weights <- function(dd, xy) {
#   # turn coordinates into SpatialPoints object:
#   # with the appropriate coordinate reference system (CRS)
#  
#    xy1 <- xy#as.matrix(xy, ncol =2, nrow = 1)#st_coordinates(xy)
#    spp <- SpatialPoints(data.frame(x = xy1[,1], y = xy1[,2]), proj4string = fm_CRS(distwb))%>%
#      st_as_sf()
#   int <- 1#st_intersection(spp, dd)
#  
#   # ddDf <- as.data.frame(dd)%>%
#   #   dplyr::mutate(X = st_coordinates(geometry)[,1],
#   #                 Y = st_coordinates(geometry)[,2])
#   weight <- log(int$weight)
#   
#   #if (any(is.na(weight))) {
#   weight[is.na(weight)] <- 0#inlabru:::bru_fill_missing(prec, spp, v$prec01)
#   #}
#   
#   
#   return(weight)
# }

weight1 <- dataBeforeClassification[[2]]%>%
  st_as_sf()
coords <- st_coordinates(weight1)
weight1New <- data.frame(x = coords[,1],
                           y = coords[,2],
                           value = weight1$weight)

weight1NewAsSf <- weight1New%>%
  st_as_sf(., coords = c("x", "y"), crs = CRS("+proj=robin +datum=WGS84 +units=km"))

ext <- c(min(mesh$loc[,1]), max(mesh$loc[,1]), min(mesh$loc[,2]), max(mesh$loc[,2]))

ras_dom <-raster(xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4],
                 crs="+proj=robin +datum=WGS84 +units=km",
                 resolution=c(3,3), vals=1)
#ras_dom
weight1 <- rasterize(weight1NewAsSf, ras_dom , "value", update = TRUE) 
weight1.spix <- as(weight1,"SpatialPixelsDataFrame")
#timeSpent.spix$layer <- timeSpent.spix$layer/1000# + 0.0001)# - mean(timeSpent.spix$layer))/sd(timeSpent.spix$layer)
#timeSpent.spix$layer <- (timeSpent.spix$layer - mean(timeSpent.spix$layer))/sd(timeSpent.spix$layer)
f.weights <- function(xy){
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  xy1 <- xy#st_coordinates(xy)
  spp <- SpatialPoints(data.frame(x = xy1[,1], y = xy1[,2]), proj4string = fm_CRS(dist))
  #proj4string(spp) <- fm_CRS(prec)
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, weight1.spix)
  #if (any(is.na(v$layer))) {
  v$layer[is.na(v$layer)] <- 1#inlabru:::bru_fill_missing(prec, spp, v$prec01)
  #}
  return(v$layer)
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




cmp1 <- list()


for(i in 1:nspecies){
  # cmp1[[i]] <- (paste0("+ betaPO0",i,
  #                      "(1)+ beta0thin(1)+cov2(main = dist, model='linear')",
  #                      "+w1",i,"(main = geometry, model =","spdes[[",i, "]])+",
  #                      "w2(main = geometry, model = spde2)+",
  #                      "cov3", i, "(main = prec, model = 'linear')",
  #                      "+cov4", i, "(main = distwb, model = 'linear')",
  #                      "+cov1",i, "(main=elev, model='linear') +beta0det",
  #                      i,"(1) +timeSpent",i, "(main=timeSpent, model='linear')",
  #                      "+ betaNina0",i,"(1)+ betaNibc", i, "(1)+betaPA", i, "(1)+ latEffect",i, "(main = decimalLatitide, model = 'iid')+yearEffect", i,"(main = year, model = 'iid')")
  # )
  
  cmp1[[i]] <- paste0("+betaNina0",i,"(1) + cov1",i, "(f.elev(st_coordinates(.data.)), model = 'linear')+cov3", i, "(f.prec(st_coordinates(.data.)), model = 'linear')+cov4", i,"(f.distwb(st_coordinates(.data.)), model = 'linear')+betaNibc0",i, "(1)+ latEffect",i, "(main = decimalLatitude, model = 'iid')+yearEffect", i,"(main = year, model = 'iid')","+beta0det",i,"(1) + timeSpent",i,"(f.time(st_coordinates(.data.)), model = 'linear')+ betaPo0",i,"(1) + w1", i, "(main = geometry, model =","spdes[[",i, "]]) + betaPA0", i, "(1) + cov2",i,"(f.dist(st_coordinates(.data.)), model = 'linear') + w2(main = geometry, model =spde2)+ weights(log(f.weights(st_coordinates(.data.))), model='offset')")
  #cmp1[[i]] <- paste0("+ betaPo0",i,"(1) + w1", i, "(main = geometry, model =","spdes[[",i, "]])")
  #cmp1[[i]] <- paste0("+betaNina0",i,"(1) + cov1",i," +cov4", i," +betaNibc0",i, "(1)+ betaPo0",i,"(1) + w1", i, "(main = geometry, model =","spdes[[",i, "]])+cov3", i)
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

spdes <- spdeslist$spdes
spde2 <- spdeslist$spde2

aa <- boundary%>%
  st_as_sf()
#ips <- ipoints(boundary, mesh)
####################################
# Fitting the model
####################################

#}

range <- diff(range(mesh$loc[,2])) 
spdeEological <- list()
#for(i in 1: nspecies){
spdeEological[[1]] <- inla.spde2.pcmatern(mesh = mesh,
                                          # PC-prior on range: P(practic.range < 0.05) = 0.01
                                          prior.range = c(30, 0.1),
                                          # PC-prior on sigma: P(sigma > 0.8) = 0.01
                                          prior.sigma = c(3, 0.1))

spdeEological[[2]] <- inla.spde2.pcmatern(mesh = mesh,
                                          # PC-prior on range: P(practic.range > range) = 0.95
                                          prior.range = c(30, 0.1),
                                          # PC-prior on sigma: P(sigma > 3) = 0.01
                                          prior.sigma = c(3, 0.1))
#}

#SPDEs for the thinning
spdeSampling <- inla.spde2.pcmatern(mesh = mesh,
                                    # PC-prior on range: P(practic.range < 0max.edge+3) = 0.01
                                    prior.range = c(30, 0.1),
                                    # PC-prior on sigma: P(sigma > 0.8) = 0.01
                                    prior.sigma = c(2, 0.1))
spdeslist <- list(spdes = spdeEological, spde2 = spdeSampling)

spdes <- spdeslist$spdes
spde2 <- spdeslist$spde2



#ips <- fm_int(mesh, samplers = aa)
# formula = as.formula(paste0("coordinates ~ betaPO0",i,"  + cov1",
#                             i,"  + cov4",
#                             i," + w1",i,"+cov3",i, "+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)+ beta0det",
#                             i,"+timeSpent",i,"+fun1(beta0det",i,", timeSpent",i,") + log(weight)")),
#Likelihoods of INLABru
likPO <- likPA <- likNina <- likNibc <- likDet <- list()
hyperNB <- list()
#Estimates of NB size from fitted model to NINA dataset
k.pm <- c(1.1197, 1.448)
for(i in 1:nspecies){
hyperNB[[i]] <- list(size = list(initial = k.pm[i],
                                 fixed = TRUE))  
}

for(i in 1:nspecies){
  likPO[[i]] <- inlabru::like("cp",
                             formula = as.formula(paste0("geometry ~ betaPo0",i,"+w1",i, "+cov1", i, "+cov3", i, "+cov4", i, "+ cov2",i,"+ weights")),
                             data = dataBeforeClassification[[i]],
                             #components = cmp,
                             domain = list(geometry = mesh),
                             #ips = fm_int(mesh, samplers = dataBeforeClassification[[i]]),
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
  # likPA[[i]] <- inlabru::like("binomial",
  #                             formula = as.formula(paste0("detdata",i," ~ betaPA0",i,"+ cov1",
  #                                                         i,"+cov3",i, "+cov4",i,"+w1",i)),
  #                             data = occurenceDf[[i]],
  #                             Ntrials = rep(1,length(detdata[[i]]$Ntrials)),
  #                             control.family= list(link = 'cloglog')#,
  #                             #ips = fm_int(mesh, samplers = detdata[[1]])#,
  #                             #components = cmp,
  #                             #domain = list(geometry = mesh),
  #                             #samplers = aa
  # )
  
  likPA[[i]] <- inlabru::like("binomial",
                              formula = as.formula(paste0("detdata",i," ~ beta0det",i," + timeSpent",i,"+betaPA0",i,"+ cov1",
                                                          i,"+cov3",i, "+cov4",i,"+w1",i)),
                              data = detdata[[i]],
                              Ntrials = detdata[[i]]$Ntrials,
                              control.family= list(link = 'cloglog'),
                              #ips = fm_int(mesh, samplers = detdata[[1]])#,
                              #components = cmp,
                              domain = list(geometry = mesh),
                              samplers = aa
  )

  likNina[[i]] <- inlabru::like("nbinomial",
                             formula = as.formula(paste0("individualCount ~ betaNina0",i, "+cov1", i, "+cov3", i, "+cov4", i,"+w1",i)),
                             data = countDf[[i]],
                             #Ntrials = detdata[[i]]$Ntrials,
                             control.family= list(hyper = hyperNB[[i]])#,
                             #ips = fm_int(mesh, samplers = countDf[[i]]),
                             #components = cmp,
                            #domain = list(geometry = mesh),
                            #samplers = aa
  )

  likNibc[[i]] <- inlabru::like("poisson",
                                formula = as.formula(paste0("individualCount ~ betaNibc0",i,"  + cov1",
                                                            i,"+cov3",i, "+cov4",i,"+w1",i)),
                                data = nbic[[i]],
                                #Ntrials = detdata[[i]]$Ntrials,
                                control.family= list(link = 'log'),
                                #ips = fm_int(mesh, samplers = nbic[[i]]),
                                #components = cmp,
                                domain = list(geometry = mesh),
                                samplers = aa
  )
}

# elev$cov11 <- elev$cov12  <- elev$alt
# dist$cov2 <- dist$distanceRaster
# prec$cov31 <- prec$cov32 <- prec$prec01
# timeSpent$timeSpent1 <- timeSpent$timeSpent2 <- timeSpent$layer
# distwb$cov41 <- distwb$cov42 <- distwb$distanceRasterWaterBody

fit5 <- inlabru::bru(cmp,
                    likPO[[1]], likPO[[2]],
                     #unquote(paste0("lik1[[", 1:nspecies, "]]", collapse = ",")),
                     #lik2[[1]],
                     #unquote(paste0("lik3[[", 1:nspecies, "]]", collapse = ",")),
                    likPA[[1]],likPA[[2]],
                  #  likDet[[1]],likDet[[2]],
                    #likNina[[1]],likNina[[2]],
                   likNibc[[1]], likNibc[[2]],
                     options = list(control.inla = list(int.strategy = "eb"#,
                                                        #control.vb=list(enable=FALSE)#,
                                            #            cmin = 0.01
                                            ),
                                    bru_method = list(
                                      #taylor = "pandemic",
                                      # search = "all",
                                      # factor = (1 + sqrt(5)) / 2,
                                      rel_tol = 0.1
                                      # max_step = 2,
                                      # lin_opt_method = "onestep"
                                    ), #change to 0.01
                                    bru_max_iter = 10,
                                    bru_verbose = 2
                                    )
)
#k.pm <- fit5$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
#k.pm2 <- fit5$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)[2]`
#kFit <- inla.emarginal(function(x) x, k.pm)
#kFit2 <- inla.emarginal(function(x) x, k.pm2)

save(fit5, file = "VSEDetectMisclassFit3.RData")

# e.int <- predict(fit5, pixels(mesh, mask = boundary), ~ exp(beta01 + cov31 + cov11+ w11 + beta0thin + cov2 + beta0det1 + timeSpent1))
# 
# ggplot() +
#   gg(e.int) +
#   gg(boundary, alpha = 0) +
#   #gg(nests, shape = "+") +
#   coord_equal()

