### Example of data generation and model fitting with Omega assumed deterministic ###
#setwd("/Users/jorgespa/Documents/Research/DataIntegration/DeadBirds/CSlibrary_new/")
## Simulating the Covariates ##


  library(spatstat)
  #Grid for simulation
  # x0 <- seq(-1.3, 4.3, length = 100)
  # y0 <- seq(-1.3,4.3, length = 100)
  # gridlocs <- expand.grid(x0,y0)
  #
  # # Covariates for true ecological state
  # gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
  # covariate.im <- im(gridcov, x0, y0)
  #
  # #Covariate for the sampling process
  # gridcov_thin <- outer(x0, y0, function(x,y) cos(2*x) - sin(2*y-4))
  # #gridcov_thin <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
  # covariate_thin.im <- im(gridcov_thin, x0, y0)
  #
  # #Covariate for the detection
  # gridcov_det <- outer(x0, y0, function(x,y) (x/2)^2+(y/2)^2)
  # covariate_detect.im <- im(gridcov_det, x0, y0)
  #
  # #Make a list for the covariates
  # cov <- list(covariate.im,covariate_thin.im,covariate_detect.im)



  ## Generating the data ##
  load("allDataForModelNew.RData")
  source("DataGeneration.R")
  library(sf)
  library(terra)
  library(raster)
  library(inlabru)
  #Extract parameters needed
  covariates = allDataForModels$covariates
  mesh = allDataForModels$mesh
  boundary = allDataForModels$boundary
  dataBeforeClassification = allDataForModels$dataBeforeClassification
  spPointsGbif = allDataForModels$spPointsGbif
  detdata = allDataForModels$detdata
  nspecies = allDataForModels$nspecies
  spdeslist = allDataForModels$spdeslist
  cov <- covariates

  elev <- covariates[[3]]
  elev$alt <- log(elev$alt)# - mean(elev$alt, na.rm = TRUE))/sd(elev$alt, na.rm = TRUE)
  dist <- covariates[[2]]
  dist$distanceRaster <- (log(dist$distanceRaster))#  - mean(log(dist$distanceRaster), na.rm = TRUE))/sd(log(dist$distanceRaster) , na.rm = TRUE)

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
                   resolution=c(3,3), vals=0.0001)
  #ras_dom
  timeSpent <- rasterize(timeSpentAsSf, ras_dom , "value", update = TRUE) 
  timeSpent$layer <- log(timeSpent)
  
  distwb <- covariates[[4]]
  distwb$distanceRasterWaterBody <- (log(distwb$distanceRasterWaterBody))# - mean(log(distwb$distanceRasterWaterBody), na.rm = TRUE))/sd(log(distwb$distanceRasterWaterBody), na.rm = TRUE)

  cov = list(elev,
             dist,
             timeSpent,
             distwb)


  nspecies <- 2#nspecies: Number of species we want to simulate
  #input: List with the parameters of the model that generates CS data
  input <-{list(
    ecological = list(
      fixed.effect=list(
        intercept = c(1, 0.5),
        betacov = c(5, 6),
        betaCov = c(-0.5, - 0.8)
      ),
      hyperparameters = list(
        sigma2 = c(2, 2),
        range = c(5, 5)
      )
    ),
    sampling = list(
      fixed.effect = list(
        intercept = c(-3),
        betacov = c(-0.2)
      ),
      hyperparameters=list(
        sigma2 = c(1),
        range = c(5)
      )
    ),
    detection = list(

      fixed.effect = list(
        intercept=c(0.8, 0.7),
        betacov = c(0.1, 0.2)
      )
    ),

    misclassification = list(

      class_prob <- matrix(c(0.9, 0.1,
                             0.05, 0.95),
                           nrow=2, ncol=2, byrow = TRUE)


    )
  )}
  # cov: a list with the covariates needed by both the data generating process and the model we want to fit.
  #idxs: A very simple list which tells which covariates belong to each stage of the data generating process
  idxs <- list(eco=c(1),
               sampling=c(2),
               detection=c(3))
  #seed: Random seed for replicating the datasets generated
  seed <- 1036610602
  #plot: Do we want to plot the results? Which results? ISSUE
  plot = list(all=FALSE,none=TRUE)
  #colmatrix: Color scale for the plots.If plot !=NULL

  # select locations
  #choose a region for use
  domain = allDataForModels$boundary
  
  ## This is a grid, useful for prediction ##
  extDom <- extent(domain)
  a <- 750#extDom@xmin
  b <- 900
  c <- 6650
  d <- 6800
  # a <- min(domain@polygons[[1]]@Polygons[[1]]@coords[,1])
  # b <- max(domain@polygons[[1]]@Polygons[[1]]@coords[,1])
  # c <- min(domain@polygons[[1]]@Polygons[[1]]@coords[,2])
  # d <- max(domain@polygons[[1]]@Polygons[[1]]@coords[,2])
  x.pred <- seq(a,b,length.out = 50)
  y.pred <- seq(c,d,length.out = 50)
  locs.pred <- expand.grid(x.pred,y.pred)
  #locs.pred <- gIntersection(domain,SpatialPoints(locs.pred,proj4string = CRS(proj4string(domain))))
  locs.pred <- SpatialPoints(locs.pred,proj4string = CRS(proj4string(domain)))
  #rast.dist <- raster("dist_hed.tif")
  
  ## Simulate a LGCP on Hedmark ##
  win <- as.owin(domain)
  x0 <- unique(locs.pred@coords[,1])
  x0 <- x0[order(x0)]
  y0 <- unique(locs.pred@coords[,2])
  y0 <- y0[order(y0)]
  
  gridlocs <- locs.pred #expand.grid(x0,y0)

library(RandomFieldsUtils)
 # seeds <- as.list(seq(1036610601, 1036610700, 1))
  seeds <- as.list(seq(1036610601, 1036610630, 1))
  simulateddata <- lapply(seeds, function(x){ csdata(nspecies=nspecies,
                          input=input,
                          cov=cov,
                          idxs=idxs,
                          seed=x,
                          plot=plot,
                          domain = allDataForModels$boundary,
                          gridlocs = gridlocs)}
  )

save(simulateddata, file = "simData.RData")

### Sampling the detections from a survey
rndpts_x0 <- runif(24, a,b)
rndpts_y0 <- runif(24, c,d)
detdata <- lapply(as.list(1:30), function(x){
  print(x)
  set.seed(1997)
  RFoptions(seed=seeds[[x]])
rndpts <- data.frame(rndpts_x0, rndpts_y0)

cov3.rast <- simulateddata[[x]]$cov[[2]] #rasterize(cov3.sp@coords,r1,cov3.sp$cov, fun=mean,na.rm=T)
cov3.spix <- as(cov3.rast,"SpatialPixelsDataFrame")

cov1.rast <- raster(simulateddata[[x]]$cov[[4]]) #rasterize(cov3.sp@coords,r1,cov3.sp$cov, fun=mean,na.rm=T)
cov1.spix <- as(cov1.rast,"SpatialPixelsDataFrame")
cov1.im <- as.im(cov1.rast)
det_prob <- occ_prob <-  list()
Eco_PP <- list()
ecoCovsIm <- list()
ecoCovsIm[[1]] <- cov1.im
## Bringing ecological covariates ##
p_eco <-1
eco_form <- c()
for(j in 1:p_eco){
  #input$ecological$fixed.effect[[",j+1,"]][i]*ecoCovsIm[[",1,"]] +
  #eco_form[j] <- paste0("input$ecological$fixed.effect[[",j+1,"]][i]*raster::extract(eco_covs, coord)")
  eco_form[j] <- paste0("input$ecological$fixed.effect[[",j+2,"]][i]*cov1.im")
}

eco_linpred <- paste(c("input$ecological$fixed.effect[[1]][i]",eco_form),collapse="+")
#eco_linpred


#win <- Window(eco_covs.im[[1]])
# Storing the raster of the true intensity for each species
#needs random field package which is currently not working
for(i in 1:nspecies){
  #x0 <- eco_covs.im[[1]]@coords[,1]
  # y0 <- eco_covs.im[[1]]@coords[,2]
  #x0 <- gridlocs@coords[,1]
  #y0 <- gridlocs@coords[,2]
  
  #x0 <- ecoCovsIm[[1]]$xcol
  #y0 <- ecoCovsIm[[1]]$yrow
  Eco_PP[[i]] <- rLGCP(model="matern",
                       eval(parse(text=eco_linpred)),
                       #lmax = 100,
                       var= input$ecological$hyperparameters[[1]][i],
                       scale= input$ecological$hyperparameters[[2]][i]/sqrt(8),
                       nu=1,
                       #win =  win,#as.owin(eval(parse(text=eco_linpred))),#owin(c(800, 850), c(6700, 6750)),
                       xy=list(x=gridlocs@coords[,1],y=gridlocs@coords[,2])
  )
  #stpp::rlgcp(s.region = win, var.grf = 1, mean.grf = 2, scale = 1)
}
eco_scale = 1

# Storing the raster of the true intensity for each species
species_rast <- list()
Scalelambda_test_rast <- list()
aa <- simulateddata[[x]]$region
for(i in 1:nspecies){
  Lam <- attr(Eco_PP[[i]], 'Lambda')
  Eco_GRF  <- log(Lam$v)
  Lambda_test <- eco_scale*Lam$v ##Test 
  # gridlocs <- gridlocs$  #expand.grid(x0,y0) ##  Matching resolutions between gridlocs and covariates
  # ecoCoords <- cbind(Eco_PP[[i]]$x, Eco_PP[[i]]$y)
  
  # df.sp <- SpatialPointsDataFrame(coords = gridlocs@coords,
  #                                 data = data.frame(w=c(log(Lam$v)),
  #                                 scalelambda =c(Lambda_test)))
  df.sp <- SpatialPointsDataFrame(coords = gridlocs,
                                  data = data.frame(w=c(Lam$v),
                                                    scalelambda =c(Lambda_test)
                                  )
  )
  
  r <- raster(df.sp)
  #r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
  xres <-  ecoCovsIm[[1]]$xstep;yres <- ecoCovsIm[[1]]$ystep## Raster resolution
  r1<-disaggregate(r, fact=res(r)/c(xres,yres))
  
  
  w1.rast <- rasterize(df.sp@coords,r1,df.sp$w, fun=mean,na.rm=T)
  w1.rastaa <- crop(w1.rast,aa)
  scalelambda.rast <- rasterize(df.sp@coords,r1,df.sp$scalelambda, fun=mean,na.rm=T) ##Test
  scalelambda.rastaa <- crop(scalelambda.rast,aa) ##Test
  
  species_rast[[i]] <- mask(w1.rastaa,aa)
  Scalelambda_test_rast[[i]] <- mask(scalelambda.rastaa,aa) ##Test
  
}
for(i in 1:nspecies){
  covs <- raster::extract(cov3.rast,rndpts)
  #covsElev <- raster::extract(Scalelambda_test_rast[[i]],rndpts)
  covsElev <- raster::extract(cov1.rast,rndpts)
  covs[is.na(covs)] <- 0
  covsElev[is.na(covsElev)] <- 0
  rndpts_lin <- rndpts %>%
    mutate(linpred = input$detection$fixed.effect$intercept[i] + input$detection$fixed.effect$betacov[i]* covs)
  rndpts_linElev <- rndpts %>%
    mutate(linpred = input$ecological$fixed.effect$intercept[i] + input$ecological$fixed.effect$betaCov[i]* covsElev)
  det_prob[[i]] <- nimble::icloglog(rndpts_lin$linpred)
  occ_prob[[i]] <- nimble::icloglog(rndpts_linElev$linpred)
}
#data_det <- rbinom(length(det_prob),1,det_prob)

#Simulate occupancy data
occData <- list()
for(i in 1:nspecies){
  occ_det <- vector("numeric", length(occ_prob[[i]]))
  for(j in 1:length(occ_prob[[i]])){
    occ_det[j]<- rbinom(1,1, occ_prob[[i]][j])
    occData[[i]] <- occ_det
  }
}


detection_data <- list()
for(i in 1:nspecies){
  data_det <- vector("numeric", length(det_prob[[i]]))
  for(j in 1:length(det_prob[[i]])){
    data_det[j]<- rbinom(1,16, occData[[i]][j]* det_prob[[i]][j])
    #data_det[j]<- rbinom(1,16,  det_prob[[i]][j])
    detection_data[[i]] <- data_det
  }
}

#Organising as spatial dataframe
data_det_spframe <- occDataSpFrame <- list()
for(i in 1:nspecies){
  filterData <- cbind(rndpts, detection_data[[i]])%>%
    dplyr::filter(`detection_data[[i]]` != 0)
  data_det_spframe[[i]] <- SpatialPointsDataFrame(cbind(filterData$rndpts_x0, filterData$rndpts_y0), data = data.frame(filterData$`detection_data[[i]]`, nTrials = 16))
  names(data_det_spframe[[i]]) <- c(paste0("detdata",i), "nTrials")
  
  occDataSpFrame[[i]] <- SpatialPointsDataFrame(rndpts,data = data.frame(occData[[i]], nTrials = 1))
  names(occDataSpFrame[[i]]) <- c(paste0("detdata",i), "nTrials")
  }
return(list(detdata1 = occDataSpFrame,
            detdata = data_det_spframe))
})

save(detdata, file = "detData.RData")

