#packages
library(ggplot2)
library(purrr)
library(dplyr)

# load truth
source("idmFunctions/functionForEstimation.R")
load("idmSimulations/simInterractionsNA50.RData")

iter = 5
allTruth <- simulations_all_na[c(1:5)]

MSEShan <- function(x, y){
  #shanIndex <- x[[2]]$summary$all.chains[ grep("shan", rownames(x[[2]]$summary$all.chains)), 1]
  #diff <- (shanIndex - log(y$hillsIDM1))
 # mahalDist <- mean(abs(diff))
 # mae <- mean(sqrt((shanIndex - y$hillsIDM1)^2))
  mae <- x[[1]]["mae",1]
  return( mae)
}


#########
# GCMSH
#######
load(paste0("idmSimulations/IG/estimateDataIGSmallshared",iter,".RData"))
gcmshShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
  purrr::flatten()%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep("GCMSH", 5),
         truth = rep("shared", 5))

gcmshPars <- lapply(rep_estimates, function(x) {
  x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
})%>%
  do.call("rbind", .) %>%
  as.data.frame()%>%
  mutate(model = rep("GCMSH", 5),
         truth = rep("shared", 5))


rm(rep_estimates)


############
# GCMCO
##########
load(paste0("idmSimulations/IG/estimateDataIGSmallcovariate",iter,".RData"))
gcmcoShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
  purrr::flatten()%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep("GCMCO", 5),
         truth = rep("shared", 5))

gcmcoPars <- lapply(rep_estimates, function(x) {
  x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
})%>%
  do.call("rbind", .) %>%
  as.data.frame()%>%
  mutate(model = rep("GCMCO", 5),
         truth = rep("shared", 5))

rm(rep_estimates)


############
# SOM
##########
load(paste0("idmSimulations/Species/estimateDataSpeciessharedSimShared",iter,".RData"))

somShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
  purrr::flatten()%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep("SOM", 5),
         truth = rep("shared", 5))

somPars <- lapply(rep_estimates, function(x) {
  x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
})%>%
  do.call("rbind", .) %>%
  as.data.frame()%>%
  mutate(model = rep("SOM", 5),
         truth = rep("shared", 5))

rm(rep_estimates)


#########
# IDMSH
#######

load(paste0("idmSimulations/IDM/estimateDataIDMsharedSimShared",iter,".RData"))

idmshShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
  purrr::flatten()%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep("IDMSH", 5),
         truth = rep("shared", 5))

idmshPars <- lapply(rep_estimates, function(x) {
  x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
})%>%
  do.call("rbind", .) %>%
  as.data.frame()%>%
  mutate(model = rep("IDMSH", 5),
         truth = rep("shared", 5))


rm(rep_estimates)


############
# IDMCO
##########

 load(paste0("idmSimulations/IDM/estimateDataIDMcovariateSimShared",iter,".RData"))

 idmcoShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
   purrr::flatten()%>%
   do.call("rbind", .)%>%
   as.data.frame()%>%
   mutate(model = rep("IDMCO", 5),
          truth = rep("shared", 5))

 idmcoPars <- lapply(rep_estimates, function(x) {
   x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
 })%>%
   do.call("rbind", .) %>%
   as.data.frame()%>%
   mutate(model = rep("IDMCO", 5),
          truth = rep("shared", 5))

 rm(rep_estimates)


#Put results together
allResSharedShan <- rbind(gcmshShan,
                          gcmcoShan,
                          somShan,
                          idmshShan,
                          idmcoShan)

allResSharedPars <- rbind(gcmshPars,
                          gcmcoPars,
                          somPars,
                          idmshPars,
                          idmcoPars)




#######################
#   Covariate
#######################

allTruth <- simulations_all_na[c(101:105)]

#########
# GCMSH
#######
load(paste0("idmSimulations/IG/estimateDataIGshared",iter,".RData"))

gcmshShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
  purrr::flatten()%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep("GCMSH", 5),
         truth = rep("covariate", 5))

gcmshPars <- lapply(rep_estimates, function(x) {
  x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
})%>%
  do.call("rbind", .) %>%
  as.data.frame()%>%
  mutate(model = rep("GCMSH", 5),
         truth = rep("covariate", 5))


rm(rep_estimates)


############
# GCMCO
##########
load(paste0("idmSimulations/IG/estimateDataIGcovariate",iter,".RData"))
gcmcoShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
  purrr::flatten()%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep("GCMCO", 5),
         truth = rep("covariate", 5))

gcmcoPars <- lapply(rep_estimates, function(x) {
  x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
})%>%
  do.call("rbind", .) %>%
  as.data.frame()%>%
  mutate(model = rep("GCMCO", 5),
         truth = rep("covariate", 5))

rm(rep_estimates)


############
# SOM
##########
load(paste0("idmSimulations/Species/estimateDataSpeciessharedSimCov",iter,".RData"))

somShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
  purrr::flatten()%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep("SOM", 5),
         truth = rep("covariate", 5))

somPars <- lapply(rep_estimates, function(x) {
  x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
})%>%
  do.call("rbind", .) %>%
  as.data.frame()%>%
  mutate(model = rep("SOM", 5),
         truth = rep("covariate", 5))

rm(rep_estimates)


#########
# IDMSH
#######

load(paste0("idmSimulations/IDM/estimateDataIDMsharedSimCov",iter,".RData"))

idmshShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
  purrr::flatten()%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep("IDMSH", 5),
         truth = rep("covariate", 5))

idmshPars <- lapply(rep_estimates, function(x) {
  x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
})%>%
  do.call("rbind", .) %>%
  as.data.frame()%>%
  mutate(model = rep("IDMSH", 5),
         truth = rep("covariate", 5))


rm(rep_estimates)


############
# IDMCO
##########
load(paste0("idmSimulations/IDM/estimateDataIDMcovariateSimCov",iter,".RData"))

idmcoShan <- purrr::map2(rep_estimates, allTruth, MSEShan)%>%
  purrr::flatten()%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep("IDMCO", 5),
         truth = rep("covariate", 5))

idmcoPars <- lapply(rep_estimates, function(x) {
  x[[1]][ grep("muSpecies|muLatitude|sigmaSpecies|sigmaSites|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda", rownames(x[[1]])), 1]
})%>%
  do.call("rbind", .) %>%
  as.data.frame()%>%
  mutate(model = rep("IDMCO", 5),
         truth = rep("covariate", 5))

rm(rep_estimates)


#Put results together
allResCovariateShan <- rbind(gcmshShan,
                          gcmcoShan,
                          somShan,
                          idmshShan,
                          idmcoShan)

allResCovariatePars <- rbind(gcmshPars,
                          gcmcoPars,
                          somPars,
                          idmshPars,
                          idmcoPars)


####################
# All data together
#####################
shanEst <- rbind(allResCovariateShan,
                  allResSharedShan)

parsEst <- rbind(allResCovariatePars,
                 allResSharedPars)

#save the data
write.csv(shanEst,paste0("simEstimates/shanSimulationsEstimates",iter,".csv"), row.names = FALSE)
write.csv(parsEst,paste0("simEstimates/parsSimulationsEstimates",iter,".csv"), row.names = FALSE)





# #load data
# gcmshShan <- extractShan(rep_estimates)
#
# gcmshPars
#
#
# #save for plotting in R
# write.csv(shanEst,"Figures/shanSimulationsEstimates.csv", row.names = FALSE)
#
#
#
#
#
#
# # load truth
# load("idmSimulations/simInterractionsNA50.RData")
# allTruth <- simulations_all_na[c(11:15, 111:115)]
#
# source("CaseStudy/functionForEstimation.R")
# load("idmSimulations/Species/estimateDataSpeciessharedSimShared.RData")
# species <- rep_estimates
#
#
#
#
# #IDMSH
# idmshSH <- lapply(as.list(0:9), function(x){
#   if(x == 0){
#     load(paste0("idmSimulations/IDM/estimateDataIDMsharedSimShared.RData"))
#   }else{
#   load(paste0("idmSimulations/IDM/estimateDataIDMsharedSimShared",x,".RData"))
#   }
#   ret <- rep_estimates
#   return(ret)
# })%>%
#   purrr::flatten()
#
# #ret <- rep_estimates
# idmshCov <- lapply(as.list(0:9), function(x){
#   if(x ==0){
#     load(paste0("idmSimulations/IDM/estimateDataIDMsharedSimCov.RData"))
#   }else{
#   load(paste0("idmSimulations/IDM/estimateDataIDMsharedSimCov",x,".RData"))
#   }
#   ret <- rep_estimates
#   return(ret)
# })%>%
#   purrr::flatten()
# #idmshCov <- c(ret, idmsh2)
#
# #put all data together
# idmsh <- c(idmshSH,
#            idmshCov)
#
#
# #IDMCO
#
# idmcoSH <- lapply(as.list(0:9), function(x){
#   if(x == 0){
#     load(paste0("idmSimulations/IDM/estimateDataIDMcovariateSimShared.RData"))
#   }else{
#     load(paste0("idmSimulations/IDM/estimateDataIDMcovariateSimShared",x,".RData"))
#   }
#   ret <- rep_estimates
#   return(ret)
# })%>%
#   purrr::flatten()
#
# #ret <- rep_estimates
# idmcoCov <- lapply(as.list(0:9), function(x){
#   if(x ==0){
#     load(paste0("idmSimulations/IDM/estimateDataIDMcovariateSimCov.RData"))
#   }else{
#     load(paste0("idmSimulations/IDM/estimateDataIDMcovariateSimCov",x,".RData"))
#   }
#   ret <- rep_estimates
#   return(ret)
# })%>%
#   purrr::flatten()
#
# idmco <- c(idmcoSH,
#            idmcoCov)
#
#
#
# # #load data
# # #IG
# # load("/Volumes/kwakupa/idmProject/idmSimulations/IG/estimateDataIGcovariate0.RData")
# # igCov1 <- rep_estimates
# # load("/Volumes/kwakupa/idmProject/idmSimulations/IG/estimateDataIGSmallcovariate0.RData")
# # igCov2 <- rep_estimates
# # igCov <- c(igCov1, igCov2)
# #
# # load("/Volumes/kwakupa/idmProject/idmSimulations/IG/estimateDataIGshared0.RData")
# # igShared1 <- rep_estimates
# # load("/Volumes/kwakupa/idmProject/idmSimulations/IG/estimateDataIGSmallshared0.RData")
# # igShared2 <- rep_estimates
# # igShared <- c(igShared1, igShared2)
# # iG <- c(igShared, igCov)
# #
# #
# # # IDM
# # load("/Volumes/kwakupa/idmProject/idmSimulations/IDM/estimateDataIDMcovariate0.RData")
# # idmCov1 <- rep_estimates
# # load("/Volumes/kwakupa/idmProject/idmSimulations/IDM/estimateDataIDMSmallcovariate0.RData")
# # idmCov2 <- rep_estimates
# # idmCov <- c(idmCov1, idmCov2)
# #
# # load("/Volumes/kwakupa/idmProject/idmSimulations/IDM/estimateDataIDMshared0.RData")
# # idmShared1 <- rep_estimates
# # load("/Volumes/kwakupa/idmProject/idmSimulations/IDM/estimateDataIDMSmallshared0.RData")
# # idmShared2 <- rep_estimates
# # idmShared <- c(idmShared1, idmShared2)
# # idm <- c(idmShared, idmCov)
# #
# #
# # # Spe
# # load("/Volumes/kwakupa/idmProject/idmSimulations/Species/estimateDataSpeciescovariate0.RData")
# # speciesCov1 <- rep_estimates
# # load("/Volumes/kwakupa/idmProject/idmSimulations/Species/estimateDataSpeciesSmallcovariate0.RData")
# # speciesCov2 <- rep_estimates
# # speciesCov <- c(speciesCov1, speciesCov2)
# #
# # load("/Volumes/kwakupa/idmProject/idmSimulations/Species/estimateDataSpeciesshared0.RData")
# # speciesShared1 <- rep_estimates
# # load("/Volumes/kwakupa/idmProject/idmSimulations/Species/estimateDataSpeciesSmallshared0.RData")
# # speciesShared2 <- rep_estimates
# # speciesShared <- c(speciesShared1, speciesShared2)
# # species <- c(speciesShared, speciesCov)
# #
# # allEstimatedData <- c(idm, iG, species)
# # # simulated data
# # load("/Volumes/kwakupa/idmProject/idmSimulations/simInterractionsNA50.RData")
# # trueShared <- simulations_all_na[1:20]
# # trueCov <- simulations_all_na[101:120]
# # truth <- c(trueShared, trueCov)
# #
# # allTruth <- rep(truth, 6)
#
# #Shannon index
#
# #load data
#
# allEstimatedData <- c(idmsh, #100
#                       idmco, #100
#                       species, #50
#                       species) #50
#
#
#
# MSEShan <- function(x, y){
#   shanIndex <- x[[2]]$summary$all.chains[ grep("shan", rownames(x[[2]]$summary$all.chains)), 1]
#   mae <- mean(sqrt((shanIndex - y$hillsIDM1)^2))
#   return(mae)
# }
#
# allTruthReplicated <- rep(allTruth, 3)
#
# shanEst <- purrr::map2(allEstimatedData, allTruthReplicated, MSEShan)%>%
#   purrr::flatten()%>%
#   do.call("rbind", .)%>%
#   as.data.frame()%>%
#   mutate(model = rep(c("IDMSH", "IDMCO", "SOM"), each = 100),
#          truth = rep(rep(c("shared", "covariate"), each = 50), 3))
#
# #save for plotting in R
# write.csv(shanEst,"Figures/shanSimulationsEstimates.csv", row.names = FALSE)
#
#
#
# #MSEPars <- function(x){
#
#   pars <- lapply(allEstimatedData, function(x) {
#     x[[2]]$summary$all.chains[ grep("muSpecies|muLatitude|muVisits|sigmaSpecies|sigmaLatitude|sigmaAlphaSpecies|sigmaAlphaSites|sigmaLambda|sigmaVisits ", rownames(x[[2]]$summary$all.chains)), 1]
#   })%>%
#   do.call("rbind", .) %>%
#     as.data.frame()%>%
#     mutate(model =c(rep("IDMSH", 99), rep("IDMCO", 99), rep("SOM", 98)),
#            truth = c(rep(c(rep("shared", 49),
#                        rep("covariate", 50)), 2),
#                      c(rep("shared", 49),
#                        rep("covariate", 49) )))
#   write.csv(pars,"Figures/parSimulationsEstimates.csv", row.names = FALSE)
#   # %>%
#   # mutate(muAlphaSpecies = muAlphaSpecies - 0,
#   #        muSpecies = muSpecies -0,
#   #        sigmaAlphaSpecies = sigmaAlphaSpecies - 1,
#   #        sigmaSpecies = sigmaSpecies - 0.2,
#   #   model = rep(c("IDM", "GCM", "SOM"), each = 80),
#   #        method = rep(rep(c("shared", "covariate"), each = 40), 3),
#   #        truth = rep(rep(c("shared", "covariate"), each = 20), 6))%>%
#   #   group_by(model, method, truth)%>%
#   #   summarise(muAlphaSpecies = mean(muAlphaSpecies),
#   #             muSpecies = mean(muSpecies),
#   #             sigmaAlphaSpecies = mean(sigmaAlphaSpecies),
#   #             sigmaSpecies = mean(sigmaSpecies)
#   #             )
#
# #save()
#
