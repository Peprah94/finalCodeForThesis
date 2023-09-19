#Script that hosts functions to format the MCMCoutput
#for summaries and plots
# x in all functions refers to the MCMCoutput

#Packages needed
library(unix)
require(dplyr)
require(ggcorrplot)
require(purrr)
require(ggpubr)
#require(heatmaply)
require(reshape2)
library(ggplot2)
require(egg)
require(readr)
library(stringr)
library(ggplot2)
#library(cowplot)
library(dplyr)
library(tidyr)
library(ggmcmc)
# devtools::install_github("https://github.com/colinharrower/BRCmap")
library(BRCmap)
if(!require(devtools)) install.packages("devtools")
#unix::rlimit_as(100*10^9)

# Extract Covariance matrix from the
#MCMC output
extractCovariance <- function(x, #x is the MCMCoutput from the analysis
                               y #y is the species names
                               ){
  nSpecies <- sqrt(length(grep("Cov", rownames(x[[3]]$summary$all.chains))))
  correlation_matrix <- matrix(x[[3]]$summary$all.chains[grep("Cov", rownames(x[[3]]$summary$all.chains)), 1],
                               nrow= nSpecies,
                               ncol= nSpecies,
                               byrow = F)%>%
    cov2cor()

  rownames(correlation_matrix) <- y
  colnames(correlation_matrix) <- y
  return(correlation_matrix)
}

# Extract Shannon index from the MCMC output
extractShan <- function(x){
  shanIndex <- matrix(x[[3]]$summary$all.chains[ grep("shan", rownames(x[[3]]$summary$all.chains)), 1],
                       nrow=1,
                       byrow = F)
  return(shanIndex)
}

# Extract Shannon SD from the MCMC output
extractShanSD <- function(x){
  shanIndex <- matrix(x[[3]]$summary$all.chains[ grep("shan", rownames(x[[3]]$summary$all.chains)), 3],
                      nrow=1,
                      byrow = F)
  return(shanIndex)
}

#function to extract pars
extractPars <- function(x){
  sigma <- matrix(x[[3]]$summary$all.chains[ grep("sigma|mu|rho", rownames(x[[3]]$summary$all.chains)), 1],
                  nrow=1,
                  #ncol=2,
                  byrow = T)
  colnames(sigma) <- rownames(x[[3]]$summary$all.chains)[grep("sigma|mu|rho", rownames(x[[3]]$summary$all.chains))]
  return(sigma)
}

# extract Sd of parameters
extractParsSD <- function(x){
  sigma <- matrix(x[[3]]$summary$all.chains[ grep("sigma|mu|rho", rownames(x[[3]]$summary$all.chains)), 3],
                  nrow=1,
                  #ncol=2,
                  byrow = T)
  colnames(sigma) <- rownames(x[[3]]$summary$all.chains)[grep("sigma|mu|rho", rownames(x[[3]]$summary$all.chains))]
  return(sigma)
}

# Function to extract and estimate hills indices for q = 0,1,2
hillsIndex <-function(x, q){
  pis = matrix((x[[13]][grep("lambda", rownames(x[[13]])), 1])[-1],
               nrow=x$n.site,
               ncol=x$n.species,
               byrow = F)%>%
    proportions(., 1)

  if(q != 1){
    hill <- (rowSums(pis^q, na.rm = TRUE))^(1/(1-q))
  }else{
    hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))
  }

  return(hill)
}

#extract Rhat
extractRhat <- function(x){
  mcmclist <- ggs(x[[3]]$samples)
  rHat <- ggs_Rhat(mcmclist, plot = FALSE)
  nRows <- nrow(rHat)
  ret <- cbind(rHat, nRows = nRows)
  return(ret)
}
