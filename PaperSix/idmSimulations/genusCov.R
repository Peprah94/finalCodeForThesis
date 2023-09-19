#Runs the analyses of the simulated data
#takes input the simulated data, methods ("IDM", "Spe", "IG"), and
#covariance_prior ("full" for sigma from invWishart), and "LV" for multiplicative
#shrinkage prior (a form of latent variable approach for sigma)
#and outputs the summary of the alpha's, beta's, z's, lambda's and correlation matrix
#as well as the values of alpha's and beta's that have converged

#
for(iCount in 1:2){
  library(doParallel)
  library(foreach)
  library(parallel)
  modelVars <- c("shared", "covariate")
  load("idmSimulations/simInterractionsNA50.RData")
  source("idmSimulations/nimbleSimulations.R")
  #source("all_sim_new.R")

  sim <- simulations_all_na[26:30]

  cl <- parallel::makeCluster(5)
  doParallel::registerDoParallel(cl)
  setDefaultCluster(cl)

  clusterExport(cl, c("run_nimble_model",
                      "incidence",
                      "nimble_incidence",
                      "richness",
                      "hill_index",
                      "nimble_hill_index",
                      "mysum",
                      "formatMatrix",
                      "nimble_sum",
                      "nimbleFormatMatrix"))

  rep_estimates <- foreach(iter = seq_along(sim), .packages = c("pbapply", "nimble", "MCMCglmm", "coda", "MCMCpack", "boot", "MASS", "parallel","ggmcmc")) %dopar% {
    tryCatch({ run_nimble_model(sim[[iter]], model = "IG",covariance_prior = "LV", method = modelVars[iCount]) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }

  save(rep_estimates, file=paste0("idmSimulations/IG/estimateDataIGSmall",modelVars[iCount],"5.RData"))
  stopCluster(cl)
}
