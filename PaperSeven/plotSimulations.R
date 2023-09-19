# Load packages
#load("gulldata/allDataForModelNew.RData")
#options(ggplot2.continuous.colour="viridis")
#options(ggplot2.continuous.fill = "viridis")

#source("DataGeneration.R")
library(sf)
library(terra)
library(raster)
library(inlabru)
library(dplyr)
library(INLA)
library(tidyr)
library(readr)
library(reshape2)
library(ggplot2)

nspecies <- 2

input <-{list(
  ecological = list(
    fixed.effect=list(
      intercept = c(-6, -5),
      betacov = c(-3, -2)
    ),
    hyperparameters = list(
      sigma2 = c(0.7, 1),
      range = c(5, 10)
    )
  ),
  sampling = list(
    fixed.effect = list(
      intercept = c(-1),
      betacov = c(1)
    ),
    hyperparameters=list(
      sigma2 = c(2),
      range = c(10)
    )
  ),
  detection = list(
    
    fixed.effect = list(
      intercept=c(0.1, 0.5),
      betacov = c(0.05, 0.1)
    )
  ),
  
  misclassification = list(
    
    class_prob <- matrix(c(0.9, 0.1,
                           0.05, 0.95),
                         nrow=2, ncol=2, byrow = TRUE)
    
    
  )
)}
#true values
input1 <-{list(
  ecological = list(
    fixed.effect=list(
      intercept = c(-5, -4),
      betacov = c(- 2, -1)
    ),
    hyperparameters = list(
      sigma2 = c(0.7, 1),
      range = c(5, 3)
    )
  ),
  sampling = list(
    fixed.effect = list(
      intercept = c( -4),
      betacov = c(-2)
    ),
    hyperparameters=list(
      sigma2 = c(1),
      range = c(5)
    )
  ),
  detection = list(
    
    fixed.effect = list(
      intercept=c(1,0.5),
      betacov = c(0.5, 0.8)
    )
  ),
  
  misclassification = list(
    
    class_prob <- matrix(c(0.9, 0.1,
                           0.05, 0.95),
                         nrow=2, ncol=2, byrow = TRUE)
    
    
  )
)}

#load data
allResults <- list()
#listsToExtract <- as.list(c(1:19, 21:30))
listsToExtract <- as.list(c(1:100))
nSamples <- length(listsToExtract)
M <- 10000
allResultsPutTogether <- pbapply::pblapply(listsToExtract, function(i){
  print(i)
load(paste0("2species example/naive/naive", i,".RData"))

  naive <- inla.posterior.sample(M, rr[[1]])
 retNaive <- lapply(naive, function(x){
   beta0Bias <- beta0RMSE <- beta1Bias <- beta1RMSE <- ecoRangeBias <- ecoSDBias <- ecoRangeRMSE <- ecoSDRMSE<- vector("numeric", 2)
    for(i in 1:nspecies){
      beta0Bias[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input$ecological$fixed.effect$intercept[i])
      beta1Bias[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input$ecological$fixed.effect$betacov[i])
      beta0RMSE[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input$ecological$fixed.effect$intercept[i])^2
      beta1RMSE[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input$ecological$fixed.effect$betacov[i])^2
      ecoRangeBias[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input$ecological$hyperparameters$range[i])
      ecoRangeRMSE[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input$ecological$hyperparameters$range[i])^2
      ecoSDBias[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input$ecological$hyperparameters$sigma2[i]))
      ecoSDRMSE[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input$ecological$hyperparameters$sigma2[i]))^2
      }
    
    gamma0Bias <- gamma0RMSE <- gamma1Bias <- gamma1RMSE <- c(NA, NA)
    alpha0Bias <- alpha0RMSE <- alpha1RMSE <- alpha1Bias <- NA
    rangeW2Bias <- rangeW2RMSE <- stdW2RMSE <- stdW2Bias <- NA

    ret <- data.frame(beta0_Bias_1 = beta0Bias[1],
                      beta0_Bias_2 = beta0Bias[2],
                      beta0_RMSE_1 = beta0RMSE[1],
                      beta0_RMSE_2 = beta0RMSE[2],
                      beta1_Bias_1 = beta1Bias[1],
                      beta1_Bias_2 = beta1Bias[2],
                      beta1_RMSE_1 = beta1RMSE[1],
                      beta1_RMSE_2 = beta1RMSE[2],
                      ecoRange_Bias_1 = ecoRangeBias[1],
                      ecoRange_Bias_2 = ecoRangeBias[2],
                      ecoRange_RMSE_1 = ecoRangeRMSE[1],
                      ecoRange_RMSE_2 = ecoRangeRMSE[2],
                      ecoSD_Bias_1 = ecoSDBias[1],
                      ecoSD_Bias_2 = ecoSDBias[2],
                      ecoSD_RMSE_1 = ecoSDRMSE[1],
                      ecoSD_RMSE_2 = ecoSDRMSE[2],
                      alpha0_Bias_1 = alpha0Bias,
                      alpha0_RMSE_1 = alpha0RMSE,
                      alpha1_Bias_1 = alpha1Bias,
                      alpha1_RMSE_1 = alpha1RMSE,
                      rangeW2_Bias_1 = rangeW2Bias,
                      rangeW2_RMSE_1 = rangeW2RMSE,
                      stdW2_Bias_1 = stdW2Bias,
                      stdW2_RMSE_1 = stdW2RMSE,
                      alpha0_Bias_2 = alpha0Bias,
                      alpha0_RMSE_2 = alpha0RMSE,
                      alpha1_Bias_2 = alpha1Bias,
                      alpha1_RMSE_2 = alpha1RMSE,
                      rangeW2_Bias_2 = rangeW2Bias,
                      rangeW2_RMSE_2 = rangeW2RMSE,
                      stdW2_Bias_2 = stdW2Bias,
                      stdW2_RMSE_2 = stdW2RMSE,
                      gamma0_Bias_1 = gamma0Bias[1],
                      gamma0_Bias_2 = gamma0Bias[2],
                      gamma0_RMSE_1 = gamma0RMSE[1],
                      gamma0_RMSE_2 = gamma0RMSE[2],
                      gamma1_Bias_1 = gamma1Bias[1],
                      gamma1_Bias_2 = gamma1Bias[2],
                      gamma1_RMSE_1 = gamma1RMSE[1],
                      gamma1_RMSE_2 = gamma1RMSE[2])
  return(ret)  
  })%>%
   do.call("rbind", .)%>%
   summarise_all(., mean, na.rm = TRUE)
  
  load(paste0("2species example/VSE/VSE", i,".RData"))
  VSE <- inla.posterior.sample(M, rr[[1]])
  retVSE <- lapply(VSE, function(x){
    beta0Bias <- beta0RMSE <- beta1Bias <- beta1RMSE <- ecoRangeBias <- ecoSDBias <- ecoRangeRMSE <- ecoSDRMSE<- vector("numeric", 2)
    gamma0Bias <- gamma0RMSE <- gamma1Bias <- gamma1RMSE <- c(NA, NA)
    alpha0Bias <- alpha0RMSE <- alpha1RMSE <- alpha1Bias <- NA
    rangeW2Bias <- rangeW2RMSE <- stdW2RMSE <- stdW2Bias <- NA
    for(i in 1:nspecies){
      beta0Bias[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input$ecological$fixed.effect$intercept[i])
      beta1Bias[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input$ecological$fixed.effect$betacov[i])
      beta0RMSE[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input$ecological$fixed.effect$intercept[i])^2
      beta1RMSE[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input$ecological$fixed.effect$betacov[i])^2
      ecoRangeBias[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input$ecological$hyperparameters$range[i])
      ecoRangeRMSE[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input$ecological$hyperparameters$range[i])^2
      ecoSDBias[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input$ecological$hyperparameters$sigma2[i]))
      ecoSDRMSE[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input$ecological$hyperparameters$sigma2[i]))^2
    }
    alpha0Bias <- (x$latent[grepl(paste0("beta0thin"), rownames(x$latent)),] - input$sampling$fixed.effect$intercept)
    alpha0RMSE <- (x$latent[grepl(paste0("beta0thin"), rownames(x$latent)),] - input$sampling$fixed.effect$intercept)^2
    alpha1RMSE <- (x$latent[grepl(paste0("cov2"), rownames(x$latent)),] - input$sampling$fixed.effect$betacov)^2
    alpha1Bias <- (x$latent[grepl(paste0("cov2"), rownames(x$latent)),] - input$sampling$fixed.effect$betacov)
    rangeW2Bias <- (x$hyperpar[grepl(paste0("Range for w2"), names(x$hyperpar))] - input$sampling$hyperparameters$range)
    rangeW2RMSE <- (x$hyperpar[grepl(paste0("Range for w2"), names(x$hyperpar))] - input$sampling$hyperparameters$range)^2
    stdW2RMSE <- (x$hyperpar[grepl(paste0("Stdev for w2"), names(x$hyperpar))] - sqrt(input$sampling$hyperparameters$sigma2))^2
    stdW2Bias <- (x$hyperpar[grepl(paste0("Stdev for w2"), names(x$hyperpar))] - sqrt(input$sampling$hyperparameters$sigma2))
    
    ret <- data.frame(beta0_Bias_1 = beta0Bias[1],
                      beta0_Bias_2 = beta0Bias[2],
                      beta0_RMSE_1 = beta0RMSE[1],
                      beta0_RMSE_2 = beta0RMSE[2],
                      beta1_Bias_1 = beta1Bias[1],
                      beta1_Bias_2 = beta1Bias[2],
                      beta1_RMSE_1 = beta1RMSE[1],
                      beta1_RMSE_2 = beta1RMSE[2],
                      ecoRange_Bias_1 = ecoRangeBias[1],
                      ecoRange_Bias_2 = ecoRangeBias[2],
                      ecoRange_RMSE_1 = ecoRangeRMSE[1],
                      ecoRange_RMSE_2 = ecoRangeRMSE[2],
                      ecoSD_Bias_1 = ecoSDBias[1],
                      ecoSD_Bias_2 = ecoSDBias[2],
                      ecoSD_RMSE_1 = ecoSDRMSE[1],
                      ecoSD_RMSE_2 = ecoSDRMSE[2],
                      alpha0_Bias_1 = alpha0Bias,
                      alpha0_RMSE_1 = alpha0RMSE,
                      alpha1_Bias_1 = alpha1Bias,
                      alpha1_RMSE_1 = alpha1RMSE,
                      rangeW2_Bias_1 = rangeW2Bias,
                      rangeW2_RMSE_1 = rangeW2RMSE,
                      stdW2_Bias_1 = stdW2Bias,
                      stdW2_RMSE_1 = stdW2RMSE,
                      alpha0_Bias_2 = alpha0Bias,
                      alpha0_RMSE_2 = alpha0RMSE,
                      alpha1_Bias_2 = alpha1Bias,
                      alpha1_RMSE_2 = alpha1RMSE,
                      rangeW2_Bias_2 = rangeW2Bias,
                      rangeW2_RMSE_2 = rangeW2RMSE,
                      stdW2_Bias_2 = stdW2Bias,
                      stdW2_RMSE_2 = stdW2RMSE,
                      gamma0_Bias_1 = gamma0Bias[1],
                      gamma0_Bias_2 = gamma0Bias[2],
                      gamma0_RMSE_1 = gamma0RMSE[1],
                      gamma0_RMSE_2 = gamma0RMSE[2],
                      gamma1_Bias_1 = gamma1Bias[1],
                      gamma1_Bias_2 = gamma1Bias[2],
                      gamma1_RMSE_1 = gamma1RMSE[1],
                      gamma1_RMSE_2 = gamma1RMSE[2])
    return(ret)  
  })%>%
    do.call("rbind", .)%>%
    summarise_all(., mean, na.rm = TRUE)
  
  
  load(paste0("2species example/VSEDetect/VSEDetect", i,".RData"))
  VSEDetect <-  inla.posterior.sample(M, rr[[1]])
  retVSEDetect <- lapply(VSEDetect, function(x){
    beta0Bias <- beta0RMSE <- beta1Bias <- beta1RMSE <- ecoRangeBias <- ecoSDBias <- ecoRangeRMSE <- ecoSDRMSE<- vector("numeric", 2)
    gamma0Bias <- gamma0RMSE <- gamma1Bias <- gamma1RMSE <- vector("numeric", 2)
    alpha0Bias <- alpha0RMSE <- alpha1RMSE <- alpha1Bias <- NA
    rangeW2Bias <- rangeW2RMSE <- stdW2RMSE <- stdW2Bias <- NA
    for(i in 1:nspecies){
      beta0Bias[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input$ecological$fixed.effect$intercept[i])
      beta1Bias[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input$ecological$fixed.effect$betacov[i])
      beta0RMSE[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input$ecological$fixed.effect$intercept[i])^2
      beta1RMSE[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input$ecological$fixed.effect$betacov[i])^2
      ecoRangeBias[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input$ecological$hyperparameters$range[i])
      ecoRangeRMSE[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input$ecological$hyperparameters$range[i])^2
      ecoSDBias[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input$ecological$hyperparameters$sigma2[i]))
      ecoSDRMSE[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input$ecological$hyperparameters$sigma2[i]))^2
      gamma0Bias[i] <- (x$latent[grepl(paste0("beta0det", i), rownames(x$latent)),] - input$detection$fixed.effect$intercept[i])
      gamma1Bias[i] <- (x$latent[grepl(paste0("cov3", i), rownames(x$latent)),] - input$detection$fixed.effect$betacov[i])
      gamma0RMSE[i] <- (x$latent[grepl(paste0("beta0det", i), rownames(x$latent)),] - input$detection$fixed.effect$intercept[i])^2
      gamma1RMSE[i] <- (x$latent[grepl(paste0("cov3", i), rownames(x$latent)),] - input$detection$fixed.effect$betacov[i])^2
    }
    alpha0Bias <- (x$latent[grepl(paste0("beta0thin"), rownames(x$latent)),] - input$sampling$fixed.effect$intercept)
    alpha0RMSE <- (x$latent[grepl(paste0("beta0thin"), rownames(x$latent)),] - input$sampling$fixed.effect$intercept)^2
    alpha1RMSE <- (x$latent[grepl(paste0("cov2"), rownames(x$latent)),] - input$sampling$fixed.effect$betacov)^2
    alpha1Bias <- (x$latent[grepl(paste0("cov2"), rownames(x$latent)),] - input$sampling$fixed.effect$betacov)
    rangeW2Bias <- (x$hyperpar[grepl(paste0("Range for w2"), names(x$hyperpar))] - input$sampling$hyperparameters$range)
    rangeW2RMSE <- (x$hyperpar[grepl(paste0("Range for w2"), names(x$hyperpar))] - input$sampling$hyperparameters$range)^2
    stdW2RMSE <- (x$hyperpar[grepl(paste0("Stdev for w2"), names(x$hyperpar))] - sqrt(input$sampling$hyperparameters$sigma2))^2
    stdW2Bias <- (x$hyperpar[grepl(paste0("Stdev for w2"), names(x$hyperpar))] - sqrt(input$sampling$hyperparameters$sigma2))
    
    ret <- data.frame(beta0_Bias_1 = beta0Bias[1],
                      beta0_Bias_2 = beta0Bias[2],
                      beta0_RMSE_1 = beta0RMSE[1],
                      beta0_RMSE_2 = beta0RMSE[2],
                      beta1_Bias_1 = beta1Bias[1],
                      beta1_Bias_2 = beta1Bias[2],
                      beta1_RMSE_1 = beta1RMSE[1],
                      beta1_RMSE_2 = beta1RMSE[2],
                      ecoRange_Bias_1 = ecoRangeBias[1],
                      ecoRange_Bias_2 = ecoRangeBias[2],
                      ecoRange_RMSE_1 = ecoRangeRMSE[1],
                      ecoRange_RMSE_2 = ecoRangeRMSE[2],
                      ecoSD_Bias_1 = ecoSDBias[1],
                      ecoSD_Bias_2 = ecoSDBias[2],
                      ecoSD_RMSE_1 = ecoSDRMSE[1],
                      ecoSD_RMSE_2 = ecoSDRMSE[2],
                      alpha0_Bias_1 = alpha0Bias,
                      alpha0_RMSE_1 = alpha0RMSE,
                      alpha1_Bias_1 = alpha1Bias,
                      alpha1_RMSE_1 = alpha1RMSE,
                      rangeW2_Bias_1 = rangeW2Bias,
                      rangeW2_RMSE_1 = rangeW2RMSE,
                      stdW2_Bias_1 = stdW2Bias,
                      stdW2_RMSE_1 = stdW2RMSE,
                      alpha0_Bias_2 = alpha0Bias,
                      alpha0_RMSE_2 = alpha0RMSE,
                      alpha1_Bias_2 = alpha1Bias,
                      alpha1_RMSE_2 = alpha1RMSE,
                      rangeW2_Bias_2 = rangeW2Bias,
                      rangeW2_RMSE_2 = rangeW2RMSE,
                      stdW2_Bias_2 = stdW2Bias,
                      stdW2_RMSE_2 = stdW2RMSE,
                      gamma0_Bias_1 = gamma0Bias[1],
                      gamma0_Bias_2 = gamma0Bias[2],
                      gamma0_RMSE_1 = gamma0RMSE[1],
                      gamma0_RMSE_2 = gamma0RMSE[2],
                      gamma1_Bias_1 = gamma1Bias[1],
                      gamma1_Bias_2 = gamma1Bias[2],
                      gamma1_RMSE_1 = gamma1RMSE[1],
                      gamma1_RMSE_2 = gamma1RMSE[2])
    return(ret)  
  })%>%
    do.call("rbind", .)%>%
    summarise_all(., mean, na.rm = TRUE)
  
  load(paste0("2species example/Detect/Detect", i,".RData"))
  Detect <- inla.posterior.sample(M, rr[[1]])
  retDetect <- lapply(Detect, function(x){
    beta0Bias <- beta0RMSE <- beta1Bias <- beta1RMSE <- ecoRangeBias <- ecoSDBias <- ecoRangeRMSE <- ecoSDRMSE<- vector("numeric", 2)
    gamma0Bias <- gamma0RMSE <- gamma1Bias <- gamma1RMSE <- vector("numeric", 2)
    alpha0Bias <- alpha0RMSE <- alpha1RMSE <- alpha1Bias <- NA
    rangeW2Bias <- rangeW2RMSE <- stdW2RMSE <- stdW2Bias <- NA
    for(i in 1:nspecies){
      beta0Bias[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input$ecological$fixed.effect$intercept[i])
      beta1Bias[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input$ecological$fixed.effect$betacov[i])
      beta0RMSE[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input$ecological$fixed.effect$intercept[i])^2
      beta1RMSE[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input$ecological$fixed.effect$betacov[i])^2
      ecoRangeBias[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input$ecological$hyperparameters$range[i])
      ecoRangeRMSE[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input$ecological$hyperparameters$range[i])^2
      ecoSDBias[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input$ecological$hyperparameters$sigma2[i]))
      ecoSDRMSE[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input$ecological$hyperparameters$sigma2[i]))^2
      gamma0Bias[i] <- (x$latent[grepl(paste0("beta0det", i), rownames(x$latent)),] - input$detection$fixed.effect$intercept[i])
      gamma1Bias[i] <- (x$latent[grepl(paste0("cov3", i), rownames(x$latent)),] - input$detection$fixed.effect$betacov[i])
      gamma0RMSE[i] <- (x$latent[grepl(paste0("beta0det", i), rownames(x$latent)),] - input$detection$fixed.effect$intercept[i])^2
      gamma1RMSE[i] <- (x$latent[grepl(paste0("cov3", i), rownames(x$latent)),] - input$detection$fixed.effect$betacov[i])^2
    }
    # alpha0Bias <- (x$latent[grepl(paste0("beta0thin"), rownames(x$latent)),] - input$sampling$fixed.effect$intercept)
    # alpha0RMSE <- (x$latent[grepl(paste0("beta0thin"), rownames(x$latent)),] - input$sampling$fixed.effect$intercept)^2
    # alpha1RMSE <- (x$latent[grepl(paste0("cov2"), rownames(x$latent)),] - input$sampling$fixed.effect$betacov)^2
    # alpha1Bias <- (x$latent[grepl(paste0("cov2"), rownames(x$latent)),] - input$sampling$fixed.effect$betacov)
    # rangeW2Bias <- (x$hyperpar[grepl(paste0("Range for w2"), names(x$hyperpar))] - input$sampling$hyperparameters$range)
    # rangeW2RMSE <- (x$hyperpar[grepl(paste0("Range for w2"), names(x$hyperpar))] - input$sampling$hyperparameters$range)^2
    # stdW2RMSE <- (x$hyperpar[grepl(paste0("Stdev for w2"), names(x$hyperpar))] - sqrt(input$sampling$hyperparameters$sigma2))^2
    # stdW2Bias <- (x$hyperpar[grepl(paste0("Stdev for w2"), names(x$hyperpar))] - sqrt(input$sampling$hyperparameters$sigma2))
    
    ret <- data.frame(beta0_Bias_1 = beta0Bias[1],
                      beta0_Bias_2 = beta0Bias[2],
                      beta0_RMSE_1 = beta0RMSE[1],
                      beta0_RMSE_2 = beta0RMSE[2],
                      beta1_Bias_1 = beta1Bias[1],
                      beta1_Bias_2 = beta1Bias[2],
                      beta1_RMSE_1 = beta1RMSE[1],
                      beta1_RMSE_2 = beta1RMSE[2],
                      ecoRange_Bias_1 = ecoRangeBias[1],
                      ecoRange_Bias_2 = ecoRangeBias[2],
                      ecoRange_RMSE_1 = ecoRangeRMSE[1],
                      ecoRange_RMSE_2 = ecoRangeRMSE[2],
                      ecoSD_Bias_1 = ecoSDBias[1],
                      ecoSD_Bias_2 = ecoSDBias[2],
                      ecoSD_RMSE_1 = ecoSDRMSE[1],
                      ecoSD_RMSE_2 = ecoSDRMSE[2],
                      alpha0_Bias_1 = alpha0Bias,
                      alpha0_RMSE_1 = alpha0RMSE,
                      alpha1_Bias_1 = alpha1Bias,
                      alpha1_RMSE_1 = alpha1RMSE,
                      rangeW2_Bias_1 = rangeW2Bias,
                      rangeW2_RMSE_1 = rangeW2RMSE,
                      stdW2_Bias_1 = stdW2Bias,
                      stdW2_RMSE_1 = stdW2RMSE,
                      alpha0_Bias_2 = alpha0Bias,
                      alpha0_RMSE_2 = alpha0RMSE,
                      alpha1_Bias_2 = alpha1Bias,
                      alpha1_RMSE_2 = alpha1RMSE,
                      rangeW2_Bias_2 = rangeW2Bias,
                      rangeW2_RMSE_2 = rangeW2RMSE,
                      stdW2_Bias_2 = stdW2Bias,
                      stdW2_RMSE_2 = stdW2RMSE,
                      gamma0_Bias_1 = gamma0Bias[1],
                      gamma0_Bias_2 = gamma0Bias[2],
                      gamma0_RMSE_1 = gamma0RMSE[1],
                      gamma0_RMSE_2 = gamma0RMSE[2],
                      gamma1_Bias_1 = gamma1Bias[1],
                      gamma1_Bias_2 = gamma1Bias[2],
                      gamma1_RMSE_1 = gamma1RMSE[1],
                      gamma1_RMSE_2 = gamma1RMSE[2])
    return(ret)  
  })%>%
    do.call("rbind", .)%>%
    summarise_all(., mean, na.rm = TRUE)
  
  load(paste0("2species example/VSEDetectMisclass/VSEDetectMisclass", i,".RData"))
  VSEDetectMisclass <- inla.posterior.sample(M, rr[[1]])
  retVSEDetectMisclass <- lapply(VSEDetectMisclass, function(x){
    beta0Bias <- beta0RMSE <- beta1Bias <- beta1RMSE <- ecoRangeBias <- ecoSDBias <- ecoRangeRMSE <- ecoSDRMSE<- vector("numeric", 2)
    gamma0Bias <- gamma0RMSE <- gamma1Bias <- gamma1RMSE <- vector("numeric", 2)
    alpha0Bias <- alpha0RMSE <- alpha1RMSE <- alpha1Bias <- NA
    rangeW2Bias <- rangeW2RMSE <- stdW2RMSE <- stdW2Bias <- NA
    for(i in 1:nspecies){
      beta0Bias[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input1$ecological$fixed.effect$intercept[i])
      beta1Bias[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input1$ecological$fixed.effect$betacov[i])
      beta0RMSE[i] <- (x$latent[grepl(paste0("beta0", i), rownames(x$latent)),] - input1$ecological$fixed.effect$intercept[i])^2
      beta1RMSE[i] <- (x$latent[grepl(paste0("cov1", i), rownames(x$latent)),] - input1$ecological$fixed.effect$betacov[i])^2
      ecoRangeBias[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input1$ecological$hyperparameters$range[i])
      ecoRangeRMSE[i] <- (x$hyperpar[grepl(paste0("Range for w1", i), names(x$hyperpar))] - input1$ecological$hyperparameters$range[i])^2
      ecoSDBias[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input1$ecological$hyperparameters$sigma2[i]))
      ecoSDRMSE[i] <- (x$hyperpar[grepl(paste0("Stdev for w1", i), names(x$hyperpar))] - sqrt(input1$ecological$hyperparameters$sigma2[i]))^2
      gamma0Bias[i] <- (x$latent[grepl(paste0("beta0det", i), rownames(x$latent)),] - input1$detection$fixed.effect$intercept[i])
      gamma1Bias[i] <- (x$latent[grepl(paste0("cov3", i), rownames(x$latent)),] - input1$detection$fixed.effect$betacov[i])
      gamma0RMSE[i] <- (x$latent[grepl(paste0("beta0det", i), rownames(x$latent)),] - input1$detection$fixed.effect$intercept[i])^2
      gamma1RMSE[i] <- (x$latent[grepl(paste0("cov3", i), rownames(x$latent)),] - input1$detection$fixed.effect$betacov[i])^2
    }
    alpha0Bias <- (x$latent[grepl(paste0("beta0thin"), rownames(x$latent)),] - input1$sampling$fixed.effect$intercept)
    alpha0RMSE <- (x$latent[grepl(paste0("beta0thin"), rownames(x$latent)),] - input1$sampling$fixed.effect$intercept)^2
    alpha1RMSE <- (x$latent[grepl(paste0("cov2"), rownames(x$latent)),] - input1$sampling$fixed.effect$betacov)^2
    alpha1Bias <- (x$latent[grepl(paste0("cov2"), rownames(x$latent)),] - input1$sampling$fixed.effect$betacov)
    rangeW2Bias <- (x$hyperpar[grepl(paste0("Range for w2"), names(x$hyperpar))] - input1$sampling$hyperparameters$range)
    rangeW2RMSE <- (x$hyperpar[grepl(paste0("Range for w2"), names(x$hyperpar))] - input1$sampling$hyperparameters$range)^2
    stdW2RMSE <- (x$hyperpar[grepl(paste0("Stdev for w2"), names(x$hyperpar))] - sqrt(input1$sampling$hyperparameters$sigma2))^2
    stdW2Bias <- (x$hyperpar[grepl(paste0("Stdev for w2"), names(x$hyperpar))] - sqrt(input1$sampling$hyperparameters$sigma2))
    
    ret <- data.frame(beta0_Bias_1 = beta0Bias[1],
                      beta0_Bias_2 = beta0Bias[2],
                      beta0_RMSE_1 = beta0RMSE[1],
                      beta0_RMSE_2 = beta0RMSE[2],
                      beta1_Bias_1 = beta1Bias[1],
                      beta1_Bias_2 = beta1Bias[2],
                      beta1_RMSE_1 = beta1RMSE[1],
                      beta1_RMSE_2 = beta1RMSE[2],
                      ecoRange_Bias_1 = ecoRangeBias[1],
                      ecoRange_Bias_2 = ecoRangeBias[2],
                      ecoRange_RMSE_1 = ecoRangeRMSE[1],
                      ecoRange_RMSE_2 = ecoRangeRMSE[2],
                      ecoSD_Bias_1 = ecoSDBias[1],
                      ecoSD_Bias_2 = ecoSDBias[2],
                      ecoSD_RMSE_1 = ecoSDRMSE[1],
                      ecoSD_RMSE_2 = ecoSDRMSE[2],
                      alpha0_Bias_1 = alpha0Bias,
                      alpha0_RMSE_1 = alpha0RMSE,
                      alpha1_Bias_1 = alpha1Bias,
                      alpha1_RMSE_1 = alpha1RMSE,
                      rangeW2_Bias_1 = rangeW2Bias,
                      rangeW2_RMSE_1 = rangeW2RMSE,
                      stdW2_Bias_1 = stdW2Bias,
                      stdW2_RMSE_1 = stdW2RMSE,
                      alpha0_Bias_2 = alpha0Bias,
                      alpha0_RMSE_2 = alpha0RMSE,
                      alpha1_Bias_2 = alpha1Bias,
                      alpha1_RMSE_2 = alpha1RMSE,
                      rangeW2_Bias_2 = rangeW2Bias,
                      rangeW2_RMSE_2 = rangeW2RMSE,
                      stdW2_Bias_2 = stdW2Bias,
                      stdW2_RMSE_2 = stdW2RMSE,
                      gamma0_Bias_1 = gamma0Bias[1],
                      gamma0_Bias_2 = gamma0Bias[2],
                      gamma0_RMSE_1 = gamma0RMSE[1],
                      gamma0_RMSE_2 = gamma0RMSE[2],
                      gamma1_Bias_1 = gamma1Bias[1],
                      gamma1_Bias_2 = gamma1Bias[2],
                      gamma1_RMSE_1 = gamma1RMSE[1],
                      gamma1_RMSE_2 = gamma1RMSE[2])
    return(ret)  
  })%>%
    do.call("rbind", .)%>%
    summarise_all(., mean, na.rm = TRUE)
 
  allResults <- rbind(retNaive,
                          retVSE,
                          retDetect,
                          retVSEDetect,
                          retVSEDetectMisclass) 
  return(allResults)
}, cl = 15)

allResults <- allResultsPutTogether%>%
  do.call("rbind", .)

samplesWithRetRes <- (dim(allResults)[1])/5

allResults <- allResults%>%
  mutate(model = rep(c("naive", "VSE", "Detect","VSEDetect", "VSEDetectMisclass"), samplesWithRetRes),
         iterations = rep(1:samplesWithRetRes, each = 5))%>%
  reshape2::melt(id.vars = c("model", "iterations"))%>%
  separate(., variable, c("Parameter", "metric", "species"))

write.csv(allResults, file = "results/simBiasRMSE.csv", row.names = FALSE)


## Plot results
simBiasRMSE <- read_csv("results/simBiasRMSE.csv")%>%
  filter(Parameter %in% c("beta1"))%>%
  mutate(value = ifelse(metric == "RMSE", sqrt(value), value))%>%
 ggplot(., mapping = aes(x = as.factor(species),
                         y = as.numeric(value),
                         fill = as.factor(model)))+
  geom_boxplot()+
  facet_wrap(~as.factor(metric), scales = "free_y")+
  theme_bw()+
  ylab("Estimates")+
  xlab("Species")+
  theme(legend.position = "bottom",
        legend.title = element_text("Model"))

#save results
ggsave("results/simBias.png",
       plot = simBiasRMSE,
       width = 7,
       height = 5,
       units = "in")

# Table of model parameters
resultsTable <- read_csv("results/simBiasRMSE.csv")%>%
  filter(!is.na(species))%>%
  mutate(value = ifelse(metric == "RMSE", sqrt(value), value))%>%
  group_by(model, metric, species, Parameter)%>%
  dplyr::select(-c(2))%>%
  summarise_all(list(mean, sd))%>%
  dplyr::rename(mean = fn1,
         sd = fn2)%>%
  reshape2::melt(id.vars = c("model", "metric",
                             "species", "Parameter"))%>%
  dplyr::rename(metricEst = variable)

#save results
write.csv(resultsTable, file = "results/simBiasRMSEtable.csv", row.names = FALSE)

## Misclassification probability
listsToExtract <- as.list(c(1:100))

misclassProbs <- lapply(listsToExtract, function(i){
load(paste0("2species example/VSEDetectMisclassProbs/VSEDetectMisclass", i,".RData"))
VSEDetectMisclass <- mcmc.out$summary[1:4,2:3]
bias <- c(VSEDetectMisclass[,1]) - c(0.9, 0.05, 0.1, 0.95)
mse <- bias <- (c(VSEDetectMisclass[,1]) - c(0.9, 0.05, 0.1, 0.95))^2
ret <- data.frame(median = c(VSEDetectMisclass[,1]),
                  sd = c(VSEDetectMisclass[,2]),
                  bias = bias,
                  iteration = i,
                  mse = mse,
                  Parameters = c("Omega_11", "Omega_21", "Omega_12", "Omega_22"))
return(ret)
})%>%
  do.call("rbind", .)%>%
  ggplot(., mapping = aes(x= Parameters, y = bias))+
  geom_boxplot()+
  theme_bw()+
  scale_x_discrete(labels = c("Omega_11" = expression(Omega[11]),
                              "Omega_12" = expression(Omega[12]),
                              "Omega_21" = expression(Omega[21]),
                              "Omega_22" = expression(Omega[22])))+
  ylim(c(0, 0.015))

#save results
ggsave("results/misclassProbs.png",
       plot = misclassProbs,
       width = 7,
       height = 5,
       units = "in")


### Extract WAIC etc

allResultsPutTogether <- pbapply::pblapply(listsToExtract, function(i){
  print(i)
  load(paste0("2species example/naive/naive", i,".RData"))
  
  naive <- c(unlist(rr[[2]][1:3]))

  
  load(paste0("2species example/VSE/VSE", i,".RData"))
  VSE <- c(unlist(rr[[2]][1:3]))

  load(paste0("2species example/VSEDetect/VSEDetect", i,".RData"))
  VSEDetect <-  c(unlist(rr[[2]][1:3]))
  
  load(paste0("2species example/Detect/Detect", i,".RData"))
  Detect <- c(unlist(rr[[2]][1:3]))
  
  load(paste0("2species example/VSEDetectMisclass/VSEDetectMisclass", i,".RData"))
  VSEDetectMisclass <- c(unlist(rr[[2]][1:3]))

  
  allResults <- rbind(naive,
                      VSE,
                      Detect,
                      VSEDetect,
                      VSEDetectMisclass) 
  return(allResults)
}, cl = 15)

allResults <- allResultsPutTogether%>%
  do.call("rbind", .)%>%
  mutate(model = rep(c("naive", "VSE", "Detect","VSEDetect", "VSEDetectMisclass"), nSamples),
         iterations = rep(1:nSamples, each = 5))%>%
  reshape2::melt(id.vars = c("model", "iterations"))%>%
  separate(., variable, c("Parameter", "metric", "species"))

write.csv(allResults, file = "results/waictable.csv", row.names = FALSE)