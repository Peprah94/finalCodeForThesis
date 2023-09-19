library(nimble)
ii <- 0
est_par <- function(p11, p22, model, return){
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
  library(terra)
  library(tidyr)
  library(stringr)
  errorIndicator <- FALSE
  listout <- list()
  omega <- matrix(c(p11, 1-p11, 1-p22, p22), 2,2, byrow = TRUE)
  if(any(omega) < 0 | any(omega) > 1) stop("Components of omega should be positive")


  csdata <- get("csdata",envir = parent.frame())
  cssampdata <- get("cssampdata",envir = parent.frame())
  detdata <- get("detdata",envir = parent.frame())
  covs <- get("covs",envir = parent.frame())
  region <- get("region",envir = parent.frame())
  #mesh <- get("mesh", envir = parent.frame())
 # listout <- get("listout",envir = parent.frame())
  ii  <- get("ii",envir =  parent.frame())
  ii <- assign("ii",ii+1,envir = parent.frame())
  crs(csdata$classifications) <- CRS("+proj=robin +datum=WGS84 +units=km")
  crs(cssampdata) <- CRS("+proj=robin +datum=WGS84 +units=km")
  

   print(ii)
  # print(x)
  # Organizing the inputs
  nspecies = nrow(omega)
  tmp <- csdata$classifications
  
  #tmp <- csdata$classification
  Eco_PPFinal_detect <- list()
  for(i in 1:nspecies){Eco_PPFinal_detect[[i]] <- tmp[which(tmp$error==i),]
  # fm_crs(Eco_PPFinal_detect[[i]]) <- fm_crs(covs[[1]])
  #slot(Eco_PPFinal_detect[[i]], "proj4string") <-  BNGproj
  crs(Eco_PPFinal_detect[[i]]) <- CRS("+proj=robin +datum=WGS84 +units=km")
  #slot(detdata[[i]], "proj4string") <-  BNGproj
  }
  # Eco_PPFinal_detect <- list()
  # for(i in 1:nspecies){
  # Eco_PPFinal_detect[[i]] <- tmp[which(tmp$error==i),]
  # slot(Eco_PPFinal_detect[[i]], "proj4string") <-  BNGproj
  # }
  Samp_PPFinal <- cssampdata
  #crs(Samp_PPFinal) <- crs(covs[[1]])
  # for(i in 1:nspecies){
  #   slot(detdata[[i]], "proj4string") <-  BNGproj
  #   #crs(detdata[[i]]) <- crs(covs[[1]])
  # }
  
  data_det_spframe <- detdata[[1]]
  detdata1 <- detdata[[2]]
  
  for(i in 1:2){
    data_det_spframe[[i]] <- data_det_spframe[[i]]
    crs(data_det_spframe[[i]]) <- CRS("+proj=robin +datum=WGS84 +units=km")
    
    detdata1[[i]] <- detdata1[[i]]
    crs(detdata1[[i]]) <- CRS("+proj=robin +datum=WGS84 +units=km")
  }

  
  #for(i in 1:2){
  #  detdata1[[i]][[1]] <- as.numeric(detdata[[i]][[1]] >0)
  #}


    #lapply(detdata, function(x){
    #crs(x) <- crs(covs[[1]])
   # return(x)
  #})

  # cov1.spix <- covs[[1]]
  # cov2.spix <- covs[[2]]
  # cov3.spix <- covs[[3]]
  aa <- region
  coords <- rbind(Eco_PPFinal_detect[[1]]@coords, Eco_PPFinal_detect[[2]]@coords)
  max.edge = diff(range(coords[,1])/(3*5))
  bound.outer <- diff(range(coords[,1])/(5))
  
  #proj4string = CRS(projection(norwaySP))
  mesh <- inla.mesh.2d(boundary = aa,
                       loc = coords,
                       max.edge  = c(1,3)*max.edge,
                       offset = c(max.edge,bound.outer),
                       cutoff = 5,
                       #min.angle = 60,
                       crs = CRS("+proj=robin +datum=WGS84 +units=km")
  )## SPDEs definition
  spdes <- list()
  for(i in 1: nspecies){
    spdes[[i]] <- inla.spde2.pcmatern(mesh = mesh,
                                      # PC-prior on range: P(practic.range < 0.05) = 0.01
                                      prior.range = c(4, 0.1),
                                      # PC-prior on sigma: P(sigma > 1) = 0.01
                                      prior.sigma = c(3, 0.1))
  }
  
  #SPDEs for the thinning
  spde2 <- inla.spde2.pcmatern(mesh = mesh,
                               # PC-prior on range: P(practic.range < 0.05) = 0.01
                               prior.range = c(4, 0.1),
                               # PC-prior on sigma: P(sigma > 1) = 0.01
                               prior.sigma = c(2, 0.1))
  
  
  # mesh <- mesh
  # spdes <- spdes$spdes
  # spde2 <- spdes$spde2


  #Defining components of the model
  cmp1 <- list()
  fun <- function(x,y,z){
    -log(1+exp((x+y+z)))
  }

  fun1 <- function(x,y){
    -log(1+exp((x+y)))
  }



  fun21 <-  function(a,b,c,d,e,f, g, h, i, j){
    ret <- log(omega[1,1]) + log(1 + ((exp(d+e+f)*plogis(i+j)*omega[2,1])/(exp(a+b+c) * plogis(g+h) *omega[1,1])))
      #log(omega[1,1]*plogis(a+b+c) + omega[2,1]*plogis(d+e+f))
    #ret <- log((exp(a+b+c)/(1+exp(a+b+c)))*omega1)
    return(ret)
  }

  fun22 <-  function(a,b,c,d,e,f, g, h, i, j){
    ret <- log(omega[2,2]) + log(1 + ((exp(a+b+c) * plogis(g+h)*omega[1,2])/(exp(d+e+f)*plogis(i+j) *omega[2,2])))
    return(ret)
  }

  # covariates
  cov1.spix$cov11 <- cov1.spix$cov12 <- cov1.spix$layer
  cov3.spix$cov31 <- cov3.spix$cov32  <- cov3.spix$layer
  crs(cov3.spix) <- crs(covs[[1]])
  names(cov2.spix)<- "cov2"
  
  f.distwb <- function(xy) {
    # turn coordinates into SpatialPoints object:
    # with the appropriate coordinate reference system (CRS)
    spp <- SpatialPoints(data.frame(x = xy[, 1], y = xy[, 2]),
                         proj4string = fm_CRS(cov1.spix)
    )
    # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
    v <- over(spp, cov1.spix)
    if (any(is.na(v$layer))) {
      v$layer <- bru_fill_missing(cov1.spix, spp, v$layer)
    }
    return(v$layer)
  }
  
  f.dist <- function(xy) {
    # turn coordinates into SpatialPoints object:
    # with the appropriate coordinate reference system (CRS)
    spp <- SpatialPoints(data.frame(x = xy[, 1], y = xy[, 2]),
                         proj4string = fm_CRS(cov2.spix)
    )
    # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
    v <- over(spp, cov2.spix)
    if (any(is.na(v$cov2))) {
      v$cov2 <- bru_fill_missing(cov2.spix, spp, v$cov2)
    }
    return(v$cov2)
  }
  
  
  f.time <- function(xy) {
    # turn coordinates into SpatialPoints object:
    # with the appropriate coordinate reference system (CRS)
    spp <- SpatialPoints(data.frame(x = xy[, 1], y = xy[, 2]),
                         proj4string = fm_CRS(cov3.spix)
    )
    # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
    v <- over(spp, cov3.spix)
    if (any(is.na(v$layer))) {
      v$cov2 <- bru_fill_missing(cov2.spix, spp, v$layer)
    }
    return(v$layer)
  }

  alpha0 <- alpha1 <- beta0 <- beta1 <- gamma0 <- gamma1<- NA

  if(model == "VSEDetect"){
  for(i in 1:nspecies){
    cmp1[[i]] <- (paste0("+ beta0",i,"(1)",
                         "+ beta0thin(1)",
                         "+ beta0det" ,i, "(1)",
                         "+w1",i,"(main = coordinates, model =","spdes[[",i, "]])",
                         "+ w2(main = coordinates, model = spde2)+",
                         "cov1",i, "(f.distwb(coordinates(.data.)),model='linear') +",
                         "cov2(f.dist(coordinates(.data.)),model='linear')+",
                         "cov3",i, "(f.time(coordinates(.data.)),model='linear')"))
  }

    cmp <- as.formula(paste0("~ -1",do.call("paste0",cmp1)))

    lik1 <- lik2 <- lik3 <- lik4 <-  list()

    for(i in 1:nspecies){
      ntrials <- detdata1[[i]]$nTrials
      lik1[[i]] <- inlabru::like("cp",
                                 formula = as.formula(paste0("coordinates ~ beta0",i,"  + cov1",i," + w1",i, "+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)+ beta0det",i,"+cov3",i,"+fun1(beta0det",i,", cov3",i,")",
                                                             "+fun2",i,"(beta01, cov11, w11, beta02, cov12, w12,beta0det1, cov31, beta0det2, cov32 )")),
                                 data = Eco_PPFinal_detect[[i]],
                                 #components = cmp,
                                 domain = list(coordinates = mesh),
                                 samplers = aa)
      lik3[[i]] <- inlabru::like("binomial",
                                 formula = as.formula(paste0("detdata",i," ~ beta0det",i,"+cov3",i)),
                                 data = detdata1[[i]],
                                 Ntrials = ntrials,#detdata1[[i]]$nTrials,
                                 control.family= list(link = 'cloglog')#,
                                 #ips = ips,
                                 #components = cmp,
                                 #domain = list(coordinates = mesh),
                                 #samplers = aa
      )
      lik2[[i]] <- inlabru::like("cp",
                                 formula = coordinates ~ beta0thin + cov2 + w2,
                                 data = Samp_PPFinal,
                                 #components = cmp,
                                 domain = list(coordinates = mesh),
                                 samplers = aa)
       lik4[[i]] <- inlabru::like("binomial",
                                  formula = as.formula(paste0("detdata",i," ~ beta0",i,"  + cov1",
                                                              i)),
                                  data = data_det_spframe[[i]],
                                  Ntrials = rep(1,24),#detdata1[[i]]$nTrials,
                                  control.family= list(link = 'cloglog')
                                  #components = cmp,
                                  #domain = list(coordinates = mesh),
                                 # samplers = aa
                                  )
    }

    #fit model
    errorIndicator <- inherits(try( fit2 <- inlabru::bru(cmp,
                         lik1[[1]], lik1[[2]],
                         lik2[[1]], #lik2[[2]],
                         lik3[[1]],lik3[[2]],
                         lik4[[1]],lik4[[2]],
                         options = list(control.inla = list(int.strategy = "eb"#,
                                                            #control.vb=list(enable=FALSE),
                                                           # cmin = 0.01
                         ),
                         #bru_method=list(rel_tol = 0.01), #change to 0.01
                         bru_max_iter =3)), 
                         silent = TRUE),
                         "try-error")

    if(errorIndicator){
      alpha0 <- gamma0 <- NA
      beta0 <- beta1 <- gamma0 <- gamma1 <- vector("numeric", 2)
      rangeW2 <- NA
      stdW2 <- NA
      rangeW11 <- NA 
      stdW11 <- NA
      rangeW12 <- NA
      stdW12 <- NA
      }else{
    tmp <- fit2$marginals.fixed
    alpha0 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("beta0thin")][[1]]))
    alpha1 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("cov2")][[1]]))
    beta0 <- beta1 <- gamma0 <- gamma1 <- vector("numeric", 2)
    for(i in 1:nspecies){
      beta0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("beta0",i)][[1]])
      beta1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]])
      gamma0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("beta0det",i)][[1]])
      gamma1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
    }

    ## The locations where information is needed
    tmpHyper <- fit2$marginals.hyperpar
    rangeW2 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w2`)
    stdW2 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w2`)
    rangeW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w11`)
    stdW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w11`)
    rangeW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w12`)
    stdW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w12`)
    }
    }

  if(model == "Detect"){
    for(i in 1:nspecies){
      cmp1[[i]] <- (paste0("+ beta0",i,"(1)",
                           "+ beta0det" ,i, "(1)",
                           "+w1",i,"(main = coordinates, model =","spdes[[",i, "]])+",
                           "cov1",i, "(f.distwb(coordinates(.data.)),model='linear') +",
                           "cov3",i, "(f.time(coordinates(.data.)),model='linear')"))
    }

    cmp <- as.formula(paste0("~ -1",do.call("paste0",cmp1)))

    lik1 <- lik2 <- lik3 <- lik4 <- list()

    for(i in 1:nspecies){
      ntrials <- detdata1[[i]]$nTrials
      lik1[[i]] <- inlabru::like("cp",
                                 formula = as.formula(paste0("coordinates ~ beta0",i,"  + cov1",i," + w1",i, "+ beta0det",i,"+cov3",i,"+fun1(beta0det",i,", cov3",i,")")),
                                 data = Eco_PPFinal_detect[[i]],
                                 #components = cmp,
                                 domain = list(coordinates = mesh),
                                 samplers = aa)
      lik3[[i]] <- inlabru::like("binomial",
                                 formula = as.formula(paste0("detdata",i," ~ beta0det",i,"+cov3",i)),
                                 data = detdata1[[i]],
                                 Ntrials = ntrials, #rep(16, length(detdata1[[i]]$detdata1)),
                                 control.family= list(link = 'cloglog')#,
                                 #ips = ips,
                                 #components = cmp,
                                 #domain = list(coordinates = mesh),
                                 #samplers = aa
      )
      lik2[[i]] <- inlabru::like("cp",
                                 formula = coordinates ~ beta0thin + cov2 + w2,
                                 data = Samp_PPFinal,
                                 #components = cmp,
                                 domain = list(coordinates = mesh),
                                 samplers = aa)
      
      lik4[[i]] <- inlabru::like("binomial",
                                 formula = as.formula(paste0("detdata",i," ~ beta0",i,"  + cov1",
                                                             i)),
                                 data = data_det_spframe[[i]],
                                 Ntrials = rep(1,24),#detdata1[[i]]$nTrials,
                                 control.family= list(link = 'cloglog')
                                 #components = cmp,
                                 #domain = list(coordinates = mesh),
                                 # samplers = aa
      )
      # lik3[[i]] <- inlabru::like("binomial",
      #                            formula = as.formula(paste0("detdata",i," ~ beta0det",i," + cov3",i)),
      #                            data = data_det_spframe[[i]],
      #                            Ntrials = data_det_spframe[[i]]$nTrials#,
      #                            #components = cmp,
      #                            #domain = list(coordinates = mesh),
      #                           # samplers = aa
      #                            )
    }

    fit2 <- inlabru::bru(cmp,
                         lik1[[1]], lik1[[2]],
                         #lik2[[1]],
                         lik3[[1]],lik3[[2]],
                         lik4[[1]],lik4[[2]],
                         options = list(control.inla = list(
                                                            int.strategy = "eb",
                                                            control.vb=list(enable=FALSE)
                                                            #cmin = 0.01
                         ),
                         #bru_method=list(rel_tol = 0.01), #change to 0.01
                         bru_max_iter =3))


    tmp <- fit2$marginals.fixed
    #alpha0 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("beta0thin")][[1]]))
    #alpha1 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("cov2")][[1]]))
    beta0 <- beta1 <- gamma0 <- gamma1 <- vector("numeric", 2)
    for(i in 1:nspecies){
      beta0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("beta0",i)][[1]])
      beta1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]])
      gamma0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("beta0det",i)][[1]])
      gamma1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]])
    }

    alpha0 <- alpha1 <- rangeW2 <- stdW2 <- NA
    ## The locations where information is needed
    tmpHyper <- fit2$marginals.hyperpar
   # rangeW2 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w2`)
    #stdW2 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w2`)
    rangeW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w11`)
    stdW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w11`)
    rangeW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w12`)
    stdW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w12`)
  }

  if(model == "VSE"){
    for(i in 1:nspecies){
      cmp1[[i]] <- (paste0("+ beta0",i,"(1)",
                           "+ beta0thin(1)",
                           "+w1",i,"(main = coordinates, model =","spdes[[",i, "]])",
                           "+ w2(main = coordinates, model = spde2)+",
                           "cov1",i, "(f.distwb(coordinates(.data.)),model='linear') +",
                           "cov2(f.dist(coordinates(.data.)),model='linear')"))
    }

    cmp <- as.formula(paste0("~ -1",do.call("paste0",cmp1)))

    lik1 <- lik2 <- lik3 <- list()

    for(i in 1:nspecies){
     # for(i in 1:nspecies){
        lik1[[i]] <- inlabru::like("cp",
                                   formula = as.formula(paste0("coordinates ~ beta0",i,"  + cov1",i," + w1",i, "+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)")),
                                   data = Eco_PPFinal_detect[[i]],
                                   #components = cmp,
                                   domain = list(coordinates = mesh),
                                   samplers = aa)

        lik2[[i]] <- inlabru::like("cp",
                                   formula = coordinates ~ beta0thin + cov2 + w2,
                                   data = Samp_PPFinal,
                                   #components = cmp,
                                   domain = list(coordinates = mesh),
                                   samplers = aa)
        # lik3[[i]] <- inlabru::like("binomial",
        #                            formula = as.formula(paste0("detdata",i," ~ beta0det",i," + cov3",i)),
        #                            data = data_det_spframe[[i]],
        #                            Ntrials = data_det_spframe[[i]]$nTrials#,
        #                            #components = cmp,
        #                            #domain = list(coordinates = mesh),
        #                           # samplers = aa
        #                            )
      }
      # lik2[[i]] <- inlabru::like("cp",
      #                            formula = coordinates ~ beta0thin + cov2 + w2,
      #                            data = Samp_PPFinal,
      #                            #components = cmp,
      #                            domain = list(coordinates = mesh),
      #                            samplers = aa)
    #}

    fit2 <- inlabru::bru(cmp,
                         lik1[[1]], lik1[[2]],
                         lik2[[1]],
                         #lik3[[1]],lik3[[2]],
                         options = list(control.inla = list(int.strategy = "eb"#,
                                                            #control.vb=list(enable=FALSE),
                                                           # cmin = 0.01
                         ),
                         #bru_method=list(rel_tol = 0.01), #change to 0.01
                         bru_max_iter =3))


    tmp <- fit2$marginals.fixed
    alpha0 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("beta0thin")][[1]]))
    alpha1 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("cov2")][[1]]))
    beta0 <- beta1 <- vector("numeric", 2)
    for(i in 1:nspecies){
      beta0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("beta0",i)][[1]])
      beta1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]])
      # gamma0 <- c(gamma0, INLA::inla.emarginal(function(x) x,tmp[paste0("beta0det",i)][[1]]))
      #gamma1 <- c(gamma1, INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]]))
    }


    gamma0 <- gamma <- c(NA, NA)
    ## The locations where information is needed
    tmpHyper <- fit2$marginals.hyperpar
    rangeW2 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w2`)
    stdW2 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w2`)
    rangeW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w11`)
    stdW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w11`)
    rangeW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w12`)
    stdW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w12`)

  }

  if(model == "naive"){
    for(i in 1:nspecies){
      cmp1[[i]] <- (paste0("+ beta0",i,"(1)",
                           "+w1",i,"(main = coordinates, model =","spdes[[",i, "]])+",
                           "cov1",i, "(f.distwb(coordinates(.data.)),model='linear') "))
    }

    cmp <- as.formula(paste0("~ -1",do.call("paste0",cmp1)))

    lik1 <- lik2 <- lik3 <- list()

    for(i in 1:nspecies){
      lik1[[i]] <- inlabru::like("cp",
                                 formula = as.formula(paste0("coordinates ~ beta0",i,"  + cov1",i," + w1",i)),
                                 data = Eco_PPFinal_detect[[i]],
                                 #components = cmp,
                                 domain = list(coordinates = mesh),
                                 samplers = aa)
      # lik3[[i]] <- inlabru::like("binomial",
      #                            formula = as.formula(paste0("detdata",i," ~ beta0",i,"  + cov1",
      #                                                        i, "beta0det",i,"+cov3",i)),
      #                            data = detdata1[[i]],
      #                            Ntrials = detdata1[[i]]$Ntrials,
      #                            control.family= list(link = 'cloglog')#,
      #                            #ips = ips,
      #                            #components = cmp,
      #                            #domain = list(coordinates = mesh),
      #                            #samplers = aa
      # )
      # lik2[[i]] <- inlabru::like("cp",
      #                            formula = coordinates ~ beta0thin + cov2 + w2,
      #                            data = Samp_PPFinal,
      #                            #components = cmp,
      #                            domain = list(coordinates = mesh),
      #                            samplers = aa)
      # lik3[[i]] <- inlabru::like("binomial",
      #                            formula = as.formula(paste0("detdata",i," ~ beta0det",i," + cov3",i)),
      #                            data = data_det_spframe[[i]],
      #                            Ntrials = data_det_spframe[[i]]$nTrials#,
      #                            #components = cmp,
      #                            #domain = list(coordinates = mesh),
      #                           # samplers = aa
      #                            )
    }

    # Fit model
    fit2 <- inlabru::bru(cmp,
                         lik1[[1]], lik1[[2]],
                         #lik2[[1]],
                         #lik3[[1]],lik3[[2]],
                         options = list(control.inla = list(
                                                            int.strategy = "eb"#,
                                                           # control.vb=list(enable=FALSE),
                                                           # cmin = 0.01
                         ),
                         #bru_method=list(rel_tol = 0.01), #change to 0.01
                         bru_max_iter =3))


    tmp <- fit2$marginals.fixed
    #alpha0 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("beta0thin")][[1]]))
    #alpha1 <- c(INLA::inla.emarginal(function(x) x,tmp[paste0("cov2")][[1]]))
    beta0 <- beta1 <- vector("numeric", 2)
    for(i in 1:nspecies){
      beta0[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("beta0",i)][[1]])
      beta1[i] <- INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]])
     # gamma0 <- c(gamma0, INLA::inla.emarginal(function(x) x,tmp[paste0("beta0det",i)][[1]]))
      #gamma1 <- c(gamma1, INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]]))
    }

    gamma0 <- gamma1 <- c(NA, NA)
    alpha0 <- alpha1 <- NA

    ## The locations where information is needed
    tmpHyper <- fit2$marginals.hyperpar
   # rangeW2 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w2`)
   # stdW2 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w2`)
    rangeW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w11`)
    stdW11 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w11`)
    rangeW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Range for w12`)
    stdW12 <- INLA::inla.emarginal(function(x) x,tmpHyper$`Stdev for w12`)
    rangeW2 = NA
    stdW2 = NA
    }



  if(errorIndicator){

  mld <- -Inf
  dic <- 0
  waic <- 0
  }else{
    mld <- fit2$mlik[1,1]
    dic <- fit2$dic$dic
    waic <- fit2$waic$waic 
}
  if(return == "inla"){
  # Get predicted values
    if(errorIndicator){
      fitVals1 <- c(rep(1, nrow(Eco_PPFinal_detect[[1]]@coords)),
                    rep(1, nrow(Eco_PPFinal_detect[[2]]@coords)))
      fitVals2 <- c(rep(1, nrow(Eco_PPFinal_detect[[1]]@coords)),
                    rep(1, nrow(Eco_PPFinal_detect[[2]]@coords)))
      
      fitTrueIntensity1 <- c(rep(1, nrow(Eco_PPFinal_detect[[1]]@coords)),
                             rep(1, nrow(Eco_PPFinal_detect[[2]]@coords)))
      fitTrueIntensity2 <- c(rep(1, nrow(Eco_PPFinal_detect[[1]]@coords)),
                             rep(1, nrow(Eco_PPFinal_detect[[2]]@coords)))
      
      
      
      #proportions
      allVals <- cbind(fitVals1,
                       fitVals2)
      props <- proportions(allVals, margin = 1)
      
      samples <- vector("numeric", length = nrow(props))
      
      for(i in 1:length(samples)){
        samples[i] <- rcat(1, props[i, ])
      }
      
      # Estimate Probability of true
      
      allVals <- cbind(fitTrueIntensity1,
                       fitTrueIntensity2)
      props <- proportions(allVals, margin = 1)
      
      fittedVals <- data.frame(lambda1 = props[,1],
                               lambda2 = props[,2],
                               samples = samples,
                               mld = -Inf,
                               dic = dic,
                               waic = waic,
                               alpha0 = alpha0,
                               alpha1 = alpha1,
                               beta01 = beta0[1],
                               beta02 = beta0[2],
                               beta11 = beta1[1],
                               beta12 = beta1[2],
                               gamma01 = gamma0[1],
                               gamma02 = gamma0[2],
                               gamma11 = gamma1[1],
                               gamma12 = gamma1[2],
                               rangeW11 = rangeW11,
                               stdW11 = stdW11,
                               rangeW12 = rangeW12,
                               stdW12 = stdW12,
                               rangeW2 = rangeW2,
                               stdW2 = stdW2)
      ret <- as.matrix(fittedVals)
    }else{
    pred <- predict(
      fit2,
      fm_pixels(mesh, mask = aa, format = "sp"),
      ~ data.frame(
        lambda1 = exp(beta01  + cov11 + w11+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)+ beta0det1+cov31+fun1(beta0det1, cov31)+fun21(beta01, cov11, w11, beta02, cov12, w12,beta0det1, cov31, beta0det2, cov32 )),
        lambda2 =exp(beta02  + cov12 + w12+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)+ beta0det2+cov32+fun1(beta0det2, cov32)+fun22(beta01, cov11, w11, beta02, cov12, w12,beta0det1, cov31, beta0det2, cov32 )),
        trueLam1 = exp(beta01  + cov11 + w11),
        trueLam2 = exp(beta02  + cov12 + w12)
      )
    )

    # Fitted values for data
    fitVals1 <- c(raster::extract(raster(pred$lambda1), Eco_PPFinal_detect[[1]]),
                  raster::extract(raster(pred$lambda1), Eco_PPFinal_detect[[2]]))
    fitVals2 <- c(raster::extract(raster(pred$lambda2), Eco_PPFinal_detect[[1]]),
                  raster::extract(raster(pred$lambda2), Eco_PPFinal_detect[[2]]))

    fitTrueIntensity1 <- c(raster::extract(raster(pred$trueLam1), Eco_PPFinal_detect[[1]]),
                  raster::extract(raster(pred$trueLam1), Eco_PPFinal_detect[[2]]))
    fitTrueIntensity2 <- c(raster::extract(raster(pred$trueLam2), Eco_PPFinal_detect[[1]]),
                  raster::extract(raster(pred$trueLam2), Eco_PPFinal_detect[[2]]))



    #proportions
    allVals <- cbind(fitVals1,
                     fitVals2)
    props <- proportions(allVals, margin = 1)

    samples <- vector("numeric", length = nrow(props))

    for(i in 1:length(samples)){
      samples[i] <- rcat(1, props[i, ])
    }

    # Estimate Probability of true

    allVals <- cbind(fitTrueIntensity1,
                     fitTrueIntensity2)
    props <- proportions(allVals, margin = 1)

    fittedVals <- data.frame(lambda1 = props[,1],
                             lambda2 = props[,2],
                             samples = samples,
                             mld = mld,
                             dic = dic,
                             waic = waic,
                             alpha0 = alpha0,
                             alpha1 = alpha1,
                             beta01 = beta0[1],
                             beta02 = beta0[2],
                             beta11 = beta1[1],
                             beta12 = beta1[2],
                             gamma01 = gamma0[1],
                             gamma02 = gamma0[2],
                             gamma11 = gamma1[1],
                             gamma12 = gamma1[2],
                             rangeW11 = rangeW11,
                             stdW11 = stdW11,
                             rangeW12 = rangeW12,
                             stdW12 = stdW12,
                             rangeW2 = rangeW2,
                             stdW2 = stdW2)
    ret <- as.matrix(fittedVals)
  }
  }else{
    fittedVals <- data.frame(  mld = mld,
                               dic = dic,
                               waic = waic,
                               alpha0 = alpha0,
                             alpha1 = alpha1,
                             beta01 = beta0[1],
                             beta02 = beta0[2],
                             beta11 = beta1[1],
                             beta12 = beta1[2],
                             gamma01 = gamma0[1],
                             gamma02 = gamma0[2],
                             gamma11 = gamma1[1],
                             gamma12 = gamma1[2],
                             rangeW11 = rangeW11,
                             stdW11 = stdW11,
                             rangeW12 = rangeW12,
                             stdW12 = stdW12,
                             rangeW2 = rangeW2,
                             stdW2 = stdW2)
    ret <- list(fit2,
                fittedVals)
  }

  return(ret)
}

#x <- 0
#a0 <- as.environment(as.list(parent.frame(), all.names=TRUE))
nimble_INLA <- nimbleRcall(
  prototype = function(
    omega=double(2) #x is a matrix
    # beta is a vector
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'est_par'
)

  #CnimbleINLA <- compileNimble(nimble_INLA)

  #Testing the compiled function.
  #Should give the same results as fit.inla

  #ii <- 0
  #listout <- list()

  class_prob <- matrix(c(0.17, 0.83,
                         0.96, 0.04),
                      nrow=2, ncol=2, byrow = TRUE)
  #CnimbleINLA(class_prob)
