load("simData.RData")
load("detData.RData")
load("allDataForModelNew.RData")
detData <- detdata
source("DataGeneration.R")
library(sf)
library(terra)
library(raster)
library(inlabru)

for(iter in 1:100){
  library(sf)
  library(terra)
  library(raster)
  library(inlabru)
  library(sp)
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


  csdata = simulateddata[[iter]]$thirdstage
  cssampdata = simulateddata[[iter]]$firststage$Samp_PPFinal
  detdata = detData[[iter]]
  covslist <- list(cov1.spix,cov2.spix,cov3.spix)
  #spdeslist <- list(spdes=spdes,spde2=spde2)
  covs = covslist
  region=simulateddata[[iter]]$region
  #mesh=simulateddata[[iter]]$mesh

  data_df <- data.frame(
    Y = csdata$classifications$error,
    C = csdata$classifications$true_species,
    eco_cov = raster::extract(cov1.rast,csdata$classifications),
    samp_cov= raster::extract(cov2.rast,csdata$classifications),
    det_cov = raster::extract(cov3.rast,csdata$classifications))


  source("estpar.R")
  tmp <- csdata$classification
  Eco_PPFinal_detect <- list()
  for(i in 1:nspecies){Eco_PPFinal_detect[[i]] <- tmp[which(tmp$error==i),]
  # fm_crs(Eco_PPFinal_detect[[i]]) <- fm_crs(covs[[1]])
  slot(Eco_PPFinal_detect[[i]], "proj4string") <-  BNGproj
  slot(detdata[[1]][[i]], "proj4string") <-  BNGproj
  slot(detdata[[2]][[i]], "proj4string") <-  BNGproj
  }
  slot(cssampdata, "proj4string") <-  BNGproj

  class_prob <- matrix(c(1, 0,
                         0, 1),
                       nrow=2, ncol=2, byrow = TRUE)
  rr <- est_par(p11 = 1, p22 = 1,
                model = "Detect",
                return = "noinla")

  save(rr, file = paste0("Detect/Detect",iter,".RData"))
}
#source("nimble_code.R")
