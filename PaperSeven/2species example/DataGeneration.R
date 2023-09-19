source("functionssimu.R")
source("Thinstage1.R")
source("Thinstage2.R")
source("Thinstage3.R")
csdata <- function(nspecies,input,cov,idxs,domain=NULL,seed,plot=list(all=TRUE),colmatrix=NULL, gridlocs){
  BNGproj <- CRS("+proj=robin +datum=WGS84 +units=km")
  #Check if the number of species provided is a number
  if(class(nspecies)!="numeric") stop("Number of species must be a number")
  
  #Setting seed for the simulation
  set.seed(seed)
  RFoptions(seed=seed)
  
  # select locations
  #choose a region for use
  # x0 <- seq(750, 900, length = 50)
  # y0 <- seq(6600, 6750 , length = 50)
  # gridlocs <- expand.grid(x0,y0)
  # coordinates(gridlocs) <- ~Var1 + Var2
  # crs(gridlocs) <- crs(domain)
  
  # The domain where the simulations are done
  # The default is provided here
  if(is.null(domain)){
    coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
    aa <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))
    win <- window(aa)#maptools::as.owin(aa) ##Need maptools
  }else{
   # domain <- crop(domain, gridlocs)
    #aa <- domain
    #win <- as.owin(aa)
    #w <- owin(c(800,850),c(6700,6750), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
    
    coordsmat <- matrix(c(a,c,b,c,b,d,a,d,a,c),ncol=2,byrow=T)
    aa <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))#%>%
      #terra::vect()
    crs(aa) <- crs(domain)
   # aaTransform <- terra::project(aa, BNGproj)%>%
    #  as(., "Spatial")

    win <- as.owin(aa)
  }
  ## Mesh for models ##
  ## Maximum distance in the extent of study area ##

  ext <- extent(aa)
  pts <- SpatialPoints(coords = matrix(ext,ncol=2),proj4string = crs(aa))
  max.dist <- spDists(pts)[1,2]

            
  
  ## Converting the covariates provided into a raster, image and pixels ##
  if(is.list(cov)){
    
    # classes <- unique(sapply(cov, class))
    # if(length(classes)==1){
    #   
    #   if(classes=="im"){
    #     covs.im <- cov
    #     covs.raster <- lapply(cov,im2rast)
    #     covs.sppixels <- lapply(covs.raster,function(x){as(x,"SpatialPixelsDataFrame")})
    #   }
    #   else{if(classes=="RasterLayer"){
    #     covs.im <- lapply(cov,as.im)
    #     covs.raster <- cov
    #     covs.sppixels <- lapply(covs.raster,function(x){as(x,"SpatialPixelsDataFrame")})
    #   }
    #     else{if(classes=="SpatialPixelsDataFrame"){
    #       covs.im <- lapply(cov,as.im)
    #       covs.raster <- lapply(cov, sppixels2raster)
    #       covs.sppixels <- cov
    #     }
    #       else{
    #         stop("Covariates must be of 'im', 'RasterLayer' or 'SpatialPixelsDataFrame'")
    #       }
    #     }}
    #   
    # }
    # else{
    #   stop("All the covariate must be in the same format")
    # }
    
    # covs.im <- lapply(cov[1:2],as.im)
    # covs.raster <- cov
    # covs.sppixels <- lapply(covs.raster,function(x){as(x,"SpatialPixelsDataFrame")})
    
    ## Indexes of the covariates##  
    eco_idxs <- idxs$eco
    sampling_idxs <- idxs$sampling
    detection_idxs <- idxs$detection
    
    eco_covs.im <- cov[eco_idxs]
    samp_covs.im <- cov[sampling_idxs]
    detect_covs.im <- cov[detection_idxs]
    
    # x <- raster(ncol=51, 
    #             nrow=51, 
    #             xmn=800, 
    #             xmx=850, 
    #             ymn=6700, 
    #             ymx=6750)
    # projection(x) <- "+proj=robin +datum=WGS84 +units=km"
    #cropElev <- crop(raster(eco_covs.im[[1]]), gridlocs)
    cropElev <- raster(eco_covs.im[[1]])
    
    #cropElev <- mask(resample(cropElev, x), aa)
    #ncol(cropElev) <- 51
    #cropElev <- mask(resample(cropDist, cropElev), aa)
    #cropElev <- aggregate(cropElev, fact=res(cropElev), FUN = mean)
    #setExtent(cropElev, ext, keepres = TRUE)
    
    #cropDist <- crop(raster(samp_covs.im[[1]]), gridlocs)
    #cropDist <- mask(resample(cropDist, cropElev), aa)
    
    cropDist <- raster(samp_covs.im[[1]])
    #cropDist <- mask(resample(cropDist, cropElev), aa)
    #res(cropDist) <- res(cropElev)
    #newres <- res(cropDist) * 1
    #cropDist <- aggregate(cropDist, fact = c(2.427799, 3.42109))
    #cropDist <- setExtent(cropDist, ext, keepres = TRUE)
    #cropTimeSpent <- crop(detect_covs.im[[1]], gridlocs)
    cropTimeSpent <- detect_covs.im[[1]]
    #setExtent(raster(cropTimeSpent), ext, keepres = TRUE)  
Distwb <- raster::raster(cov[[4]])

     ecoCovsIm <- list()
     
     cropElev <- crop(cropElev,coordsmat)
     #cropElev <- mask(cropElev,domain)
     ecoCovsIm[[1]] <- as.im(cropElev)
     #crs(cropElev)
     #Window(ecoCovsIm[[1]]) <- win
     
     cropDistwb <- crop(Distwb, coordsmat)
     cropDistwbIm <- as.im(cropDistwb)
     #cropDistwb <- mask(cropDistwb, domain)
ecoCovsIm[[2]] <- cropDistwbIm
    
    if(!is.null(input$ecological$fixed.effect$scale)){eco_scale=input$ecological$fixed.effect$scale} ##Test
    else{eco_scale=1} ##Test
    if(!is.null(input$sampling$fixed.effect$scale)){sampling_scale=input$sampling$fixed.effect$scale}
    else{sampling_scale=1}
    if(!is.null(input$detection$fixed.effect$scale)){detection_scale=input$detection$fixed.effect$scale}
    else{detection_scale=1}
    
    ## Generate ecological process ##
    
    #Ecological Process
    Eco_PP <- list()
    
    ## Bringing ecological covariates ##
    p_eco <- length(eco_idxs)
    eco_form <- c()
    for(j in 1:p_eco){
      #input$ecological$fixed.effect[[",j+1,"]][i]*ecoCovsIm[[",1,"]] +
      #eco_form[j] <- paste0("input$ecological$fixed.effect[[",j+1,"]][i]*raster::extract(eco_covs, coord)")
      eco_form[j] <- paste0("input$ecological$fixed.effect[[",j+2,"]][i]*ecoCovsIm[[",2,"]]")
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
  
    
    # Storing the raster of the true intensity for each species
    species_rast <- list()
    Scalelambda_test_rast <- list()
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

      #) ##Test
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
    
    #environment_list <- as.list(environment())
    #### First thinning stage ##
    #print(length(environment_list))
    sampling_covs.im <- list()
    cropDist <- crop(cropDist,coordsmat)
    #cropDist <- mask(cropDist,domain)
    #cropDist <- resample(raster(cropDist), cropElev)
    sampling_covs.im[[1]] <- as.im(cropDist[[1]])

    #Window(sampling_covs.im[[1]]) <- win
    firststage <- firstthinning(input,scale=sampling_scale)
    
    ## Second thinning stage ##
    detection_covs.raster <- list()
    cropTimeSpent <- crop(cropTimeSpent,coordsmat)
    #cropTimeSpent <- mask(raster(cropTimeSpent),domain)
    #ecoCovsIm[[1]] <- as.im(cropElev)
    #cropTimeSpent <- resample(raster(cropTimeSpent), Scalelambda_test_rast[[1]])
    detection_covs.raster[[1]] <- cropTimeSpent
    environment_list <- as.list(environment())
    secondstage <- secondthinning(input,environment_list,scale=detection_scale)
    
    ## Second thinning stage ##
    environment_list <- as.list(environment())
    thirdstage <- thirdthinning(input,environment_list)
    
    if(plot$all==TRUE){
      ##True ecological ##
      lapply(1:nspecies,function(x){plot(species_rast[[x]],axes=F,box=F,main=paste0("True Ecological State for species",x))
        points(Eco_PP[[x]],pch=19,cex=0.5) })
      
      ## Detection##
      
      lapply(1:nspecies,function(x){plot(secondstage$detectionprobraster[[x]],axes=F,box=F,main=paste0("Detection probability for species ",x))
        points(firststage$Eco_PPFinal[[x]],pch=19,cex=0.5)
        points(secondstage$Eco_PPFinal_detect[[x]],pch=19,cex=0.3,col="red") })
      
      ## Sampling ##
      lapply(1:nspecies,function(x){plot(firststage$retainprobraster,axes=F,box=F,main=paste0("Retaining probability for species ",x))
        points(Eco_PP[[x]],pch=19,cex=0.5)
        points(firststage$Eco_PPFinal[[x]],pch=19,cex=0.3,col="red")
      })
      
      ## Classification ##
      if(is.null(colmatrix)){
        colmatrix <- matrix(NA,nrow=nspecies,ncol=2)
        colmatrix[,1] <- sample(colors(),size = nspecies,replace = FALSE)
        colmatrix[,2] <- as.numeric(sample(0:25,size = nspecies,replace = FALSE))
        colmatrix <- data.frame(colmatrix)
        names(colmatrix) <- c("color","pch")
        colmatrix$pch <- as.numeric(colmatrix$pch)
      }
      
      for(i in 1:nspecies){
        plot(aa,axes=F,main=paste0("CS reports for species ",i))
        points(thirdstage$classifications[which(thirdstage$classifications$true_species==i),],pch=colmatrix[i,2],cex=1,col=colmatrix[i,1])
        indexes0 <- 1:nspecies
        indexes <- indexes0[-indexes0[which(indexes0==i)]]
        for(j in indexes){
          try(points(thirdstage$classifications[which(thirdstage$classifications$true_species==i & thirdstage$classifications$error==j),],pch=colmatrix[j,2],cex=0.6,col=colmatrix[j,1]))
        }
        legendtext <- paste0("Species ",1:nspecies)
        legend("right", legend=legendtext,col=colmatrix$color,pch=colmatrix$pch)
      }
      
    }
    else{
      plotnew <- within(plot,rm(all))
      which.plot <- names(plotnew[sapply(plotnew,isTRUE)])
      
      if("ecological"%in%which.plot){
        lapply(1:nspecies,function(x){plot(species_rast[[x]],axes=F,box=F,main=paste0("True Ecological State for species",x))
          points(Eco_PP[[x]],pch=19,cex=0.5)
        })
      }
      
      if("detection"%in%which.plot){
        lapply(1:nspecies,function(x){plot(secondstage$detectionprobraster[[x]],axes=F,box=F,main=paste0("Detection probability for species ",x))
          points(firststage$Eco_PPFinal[[x]],pch=19,cex=0.5)
          points(secondstage$Eco_PPFinal_detect[[x]],pch=19,cex=0.3,col="red")
          
        })
        
      }
      
      if("sampling"%in%which.plot){
        lapply(1:nspecies,function(x){plot(firststage$retainprobraster,axes=F,box=F,main=paste0("Retaining probability for species ",x))
          points(Eco_PP[[x]],pch=19,cex=0.5)
          points(firststage$Eco_PPFinal[[x]],pch=19,cex=0.3,col="red")
        })
      }
      
      if("classification"%in%which.plot){
        if(is.null(colmatrix)){
          colmatrix <- matrix(NA,nrow=nspecies,ncol=2)
          colmatrix[,1] <- sample(colors(),size = nspecies,replace = FALSE)
          colmatrix[,2] <- as.numeric(sample(0:25,size = nspecies,replace = FALSE))
          colmatrix <- data.frame(colmatrix)
          names(colmatrix) <- c("color","pch")
          colmatrix$pch <- as.numeric(colmatrix$pch)
        }
        
        for(i in 1:nspecies){
          plot(aa,axes=F,main=paste0("CS reports for species ",i))
          points(thirdstage$classifications[which(thirdstage$classifications$true_species==i),],pch=colmatrix[i,2],cex=1,col=colmatrix[i,1])
          indexes0 <- 1:nspecies
          indexes <- indexes0[-indexes0[which(indexes0==i)]]
          for(j in indexes){
            try(points(thirdstage$classifications[which(thirdstage$classifications$true_species==i & thirdstage$classifications$error==j),],pch=colmatrix[j,2],cex=0.6,col=colmatrix[j,1]))
          }
          legendtext <- paste0("Species ",1:nspecies)
          legend("right", legend=legendtext,col=colmatrix$color,pch=colmatrix$pch)
        }
        
      }
      
    }
    
  }
  else{
    stop("Covariates input should be a list")
  }
  
  mesh <- inla.mesh.2d(loc.domain = firststage$Samp_PPFinal,
                       max.edge = c(100,500),
                       offset = c(10, 50),
                       boundary = domain,
                       crs = crs(aa),
                       cutoff = 5 #0.2*max.dist
  )
  
  ## Generating lambda_obs ##
  # lambda_obs_raster <- list() ##Test
  # for(i in 1:nspecies){ ##Test
  #   retainProbs <- resample(raster(firststage$retainprobraster), Scalelambda_test_rast[[1]])
  #   detectProbs <-resample(raster(secondstage$detectionprobraster[[i]]), Scalelambda_test_rast[[1]])
  #   lambda_obs_raster[[i]] <- Scalelambda_test_rast[[i]]*(retainProbs)*(detectProbs)
  # }
  
  cov <- list(cropElev,
              cropDist[[1]], 
              cropTimeSpent,
              cropDistwbIm)
  
  return(list(trueecological=Eco_PP,firststage=firststage,secondstage=secondstage,thirdstage=thirdstage,
              species_raster=species_rast,
              #lambda_obs_raster=lambda_obs_raster,
              cov = cov,
              mesh = mesh,
              region = aa,
              input = input))
  
}

