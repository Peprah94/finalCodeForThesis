#Parameters needed for this function
load("/Volumes/kwakupa/idmProject/idmFunctions/parameters_50.RData")
source("/Volumes/kwakupa/idmProject/idmFunctions/functionForSimulationWithMissing.R")
#Simulation of the data with replication at sites

#List of input parameters
input_list_na <- list(input50)

#Number of replicates
nreplicates <- 100
simulationsShared <- pblapply(input_list_na, function(x){
  pblapply(1:nreplicates, function(z){
    data <- sim(x, seed = z, model = "shared")
  }, cl=4)
}, cl=4)

simulationsCovariate <- pblapply(input_list_na, function(x){
  pblapply(1:nreplicates, function(z){
    data <- sim(x, seed = z, model = "covariate")
  }, cl=4)
}, cl=4)



simulations_all_na <- flatten(c(simulationsShared,
                                simulationsCovariate))

#save the results
save(simulations_all_na, file="simInterractionsNA50.RData")
save(input_list_na, file="sim_input_na.RData")



