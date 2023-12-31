## Second stage of thinning (detection) ##

secondthinning <- function(input, ...){
  environment_list <- as.list(parent.frame())
  p_detection <- length(environment_list$detection_idxs)
  detection_form <- c()
  
  for(j in 1:p_detection){
    detection_form[j] <- paste0("input$detection$fixed.effect[[",j+1,"]][i]*raster::extract(environment_list$detection_covs.raster[[",j,"]],environment_list$firststage$Eco_PPFinal[[i]])")
  }  
  
  detection_linpred <- paste(c("input$detection$fixed.effect[[1]][i]",detection_form),collapse="+")
  
  #raster::extract(raster(environment_list$detection_covs.raster[[1]]),environment_list$firststage$Eco_PPFinal[[1]])
  
  # Computing detection probabilities
  det_PP <- list()
  for(i in 1:nspecies){
    det_PP[[i]] <-  psych::logistic(eval(parse(text=detection_linpred)))
  }

  #Point pattern from thinning due to detection
  Eco_PPFinal_detect <- list()
  for(i in 1:nspecies){
    if(length(environment_list$firststage$Eco_PPFinal[[i]])>0){
    environment_list$firststage$Eco_PPFinal[[i]]$retain_detect <- sapply(1:length(environment_list$firststage$Eco_PPFinal[[i]]),
                                             function(x){rbinom(1,1,p=det_PP[[i]][x])})
    Eco_PPFinal_detect[[i]] <- environment_list$firststage$Eco_PPFinal[[i]][which(environment_list$firststage$Eco_PPFinal[[i]]$retain_detect==1),]
    }
    else{Eco_PPFinal_detect[[i]] <- environment_list$firststage$Eco_PPFinal[[i]]}
  }
  
  ## Making rasters of detection probability ##

  
  detrasts <- list()
  detectionraster_form <- c()
  
  for(j in 1:p_detection){
    detectionraster_form[j] <- paste0("input$detection$fixed.effect[[",j+1,"]][i]*environment_list$detection_covs.raster[[",j,"]]")
  }  
  
  detectionraster_linpred <- paste(c("input$detection$fixed.effect[[1]][i]",detectionraster_form),collapse="+")
  
  for(i in 1:nspecies){
    detrasts[[i]] <- psych::logistic(eval(parse(text=detectionraster_linpred)))
  }
  
  detrast.aa <- lapply(detrasts,function(x){
    detrast.aa <- crop(x,environment_list$aa)
    detrast.aa <- mask(detrast.aa,environment_list$aa)
  })
  
  
  
return(list(Eco_PPFinal_detect=Eco_PPFinal_detect,
            detectionprobraster = detrast.aa))  

}

