#setwd("/Volumes/kwakupa/idmProject/CaseStudy")
load("CaseStudy/data_idm.RData")
source("idmFunctions/functionForEstimation.R")

#Extracting values from formatted data for Table 1 in main paper
lapply(simulations_all, function(x){
  sumIsNA <- function(z){
  sum(!is.na(c(z)))/nrow(z)
  }
  missSitesSpecies <- sum(x$mat.species > 0 , na.rm = TRUE)/length(!is.na(x$mat.species))
  missSitesSCounts <- sum(x$mat.genus, na.rm = TRUE)/length(!is.na(x$mat.genus))
  missSitesSpeciesSD <- sd(x$mat.species > 0, na.rm = TRUE)
  missSitesCountsSD <- sd(x$mat.genus, na.rm = TRUE)
 ret <- c(missSitesSpecies,  missSitesSCounts,
          missSitesSpeciesSD, missSitesCountsSD)
 return(ret)
   })



# load data

#DataFormatting <- function(name){

#loading data from folder
# loadData <- function(name){
#   if(!name %in% c("bumblebees", "hoverflies", "solitarybees")) stop("name must be either 'bumbebees', 'hoverflies' or 'solitarybees'.")
#   load(paste("CaseStudy",name,"estimateCaseStudyIDM.RData", sep="/"))
#   idm <- estimatesCaseStudy
#     load(paste("CaseStudy",name,"estimateCaseStudyIG.RData", sep="/"))
#     groupCnts <- estimatesCaseStudy
#   load(paste("CaseStudy",name,"estimateCaseStudySpe.RData", sep="/"))
#   speciesOcc <- estimatesCaseStudy
#
#   ret <- c(idm, groupCnts, speciesOcc)
#
#   rm(idm, groupCnts, speciesOcc, estimatesCaseStudy)
#   return(ret)
# }

insectGroups = list("bumblebees", "hoverflies", "solitarybees")
allData <-  lapply(as.list(c("bumblebees", "hoverflies", "solitarybees")), function(x){
  load(paste0("CaseStudy/",x,"/estimateCaseStudy.RData"))
  idm <- estimatesCaseStudy

  load(paste0("CaseStudy/",x,"/estimateCaseStudyIG.RData"))
  genus <- estimatesCaseStudy

  ret <- c(idm,
           genus)

   return(ret)
})

#Put all data together
allDataFlat <- allData%>%
  purrr::flatten(.)

#order of the data
# Index   |   Insect Group    |   Method    |   Type
# 1               BB                IDM       Shared
# 2                                 IDM       Covariate
# 3                                 GO        Shared
# 4                                 GO        Covariate
# 5                                 Spe       Shared
# 6                                 Spe       Covariate
# 7               HV                IDM       Shared
# 8                                 IDM       Covariate
# 9                                 GO        Shared
# 10                                GO        Covariate
# 11                                Spe       Shared
# 12                                Spe       Covariate
# 13               SB                IDM       Shared
# 14                                 IDM       Covariate
# 15                                 GO        Shared
# 16                                GO        Covariate
# 17                                Spe       Shared
# 18                                Spe       Covariate
#save(allDataFlat, file = "allDataPutTogether.RData")

# Retrieving names of species
speciesNamesList <- c(rep(list(simulations_all[[1]]$species_names),5),
                      rep(list(simulations_all[[2]]$species_names),5),
                      rep(list(simulations_all[[3]]$species_names),5)
)


# list of constants
nGroups = 3 #number of insect groups
nModels = 5 #number of models
nMethods = 3


#retrieving constants and adding it to the MCMCoutput data
simData = rep(simulations_all, each = nModels) #Data used for MCMC


#Abbreviate the names of the species
abbreSpeciesNamesList <- lapply(speciesNamesList, function(x){
  ret <- abbreviate(x, minlength = 8)
  sub("(\\w+\\s+\\w+).*", "\\1", ret)
})


##################
#   Species Interraction plot
##############
# correlationMatrices <- pmap(list(allDataFlat,
#                                  abbreSpeciesNamesList),
#                             extractCovariance)
#
# save(correlationMatrices, file = "Figures/correlationMatrices.RData")

##########################################
# extracting the Fixed effect and random effect parameters
############################################
#extract the posterior mean
parsEstimates <- lapply(allDataFlat, extractPars)%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO", "IGSH", "IGCO"), nGroups),
         insects = rep(c("Bumblebees", "Hoverflies", "Solitary bees"),
                       each= nModels))%>%
  reshape2::melt(id.vars = c("insects", "model"))

#extract the sd of the estimates
parsEstimatesSD <- lapply(allDataFlat, extractParsSD)%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO", "IGSH", "IGCO"), nGroups),
         insects = rep(c("Bumblebees", "Hoverflies", "Solitary bees"),
                       each= nModels))%>%
  reshape2::melt(id.vars = c("insects", "model"))

#extract lower CI
parsEstimatesLower <- lapply(allDataFlat, extractParsLower)%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO", "IGSH", "IGCO"), nGroups),
         insects = rep(c("Bumblebees", "Hoverflies", "Solitary bees"),
                       each= nModels))%>%
  reshape2::melt(id.vars = c("insects", "model"))

#extract upper CI
parsEstimatesUpper <- lapply(allDataFlat, extractParsUpper)%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO", "IGSH", "IGCO"), nGroups),
         insects = rep(c("Bumblebees", "Hoverflies", "Solitary bees"),
                       each= nModels))%>%
  reshape2::melt(id.vars = c("insects", "model"))

parsData <- cbind(parsEstimates,
                  sd = parsEstimatesSD$value,
                  lower = parsEstimatesLower$value,
                  upper = parsEstimatesUpper$value)

#colnames(sigmaEstimates) <- c("Site detection", "Site Observed", "Methods", "Insect Group")
write.csv(parsData,"Figures/sigmaEstimates.csv", row.names = FALSE)

####################
# Shannon Diversity
#####################

# Mean shannon Index
shannonEstimates <- lapply(allDataFlat, extractShan)%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO", "IGSH", "IGCO"), nGroups),
         insects = rep(c("Bumblebees", "Hoverflies", "Solitary bees"),
                       each= nModels))%>%
  reshape2::melt(id.vars = c("insects", "model"))%>%
  mutate(sites = rep(simulations_all[[1]]$sites$X1km_square, each = nGroups*nModels))

#colnames(shannonEstimates) = c(simulations_all[[1]]$sites$X1km_square, "insects", "model")

#shannonEstimatesFlat <- shannonEstimates%>%
#  reshape2::melt(id.vars = c( "insects", "model"))


# Standard deviation for Shannon estimates
shannonEstimatesSD <- lapply(allDataFlat, extractShanSD)%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO", "IGSH", "IGCO"), nGroups),
         insects = rep(c("Bumblebees", "Hoverflies", "Solitary bees"),
                       each= nModels))%>%
  reshape2::melt(id.vars = c("insects", "model"))%>%
  mutate(sites = rep(simulations_all[[1]]$sites$X1km_square, each = nGroups*nModels))

#colnames(shannonEstimatesSD) = c(simulations_all[[1]]$sites$X1km_square, "method", "insects", "model")

#shannonEstimatesFlatSD <- shannonEstimatesSD%>%
 # reshape2::melt(id.vars = c("method", "insects", "model"))

allShannonEst <- cbind(shannonEstimates,  sd = shannonEstimatesSD$value )
colnames(allShannonEst)[4] <- c("mean")

#add coordinates to the Shannon data
#coords <- gr_let2num(allShannonEst$sites)
coords <- allShannonEst$sites
allDataShannon <- cbind(coords, allShannonEst)

#save shannon data
write.csv(allDataShannon,"Figures/shannonEstimates.csv", row.names = FALSE)


##################
# Convergence Parameters
#################
rHatEstimatesSD <- lapply(allDataFlat, extractRhat)%>%
  do.call("rbind", .)%>%
  as.data.frame()

write.csv(rHatEstimatesSD,"Figures/rHatEstimates.csv", row.names = FALSE)


##########
# Extract WAIC
########

waicEstimates <- lapply(allDataFlat, extractWAIC)%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO", "IGSH", "IGCO" ), nGroups),
         insects = rep(c("Bumblebees", "Hoverflies", "Solitary bees"),
                       each= nModels))

write.csv(waicEstimates,"Figures/waicEstimates.csv", row.names = FALSE)
