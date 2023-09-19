# This script runs the "constant" model with all the simulated data
# Takes as input the simulated data, and run_simulations function

#load required packages
library(doParallel)

#Run for each parameter configurations
for(iCount in 1:16){
  message(paste("Running code iteration:", iCount))
  #load data
  load(paste0("simData/simulatedData",iCount,".RData"))

  #load function needed to run the MCMC in NIMBLE
  source("all_sim_new.R")

  #the data
  sim <- sim[1:600]

  # Cluster
  cl <- parallel::makeCluster(10)
  doParallel::registerDoParallel(cl)
  setDefaultCluster(cl)

  #export functions to the cluster
  clusterExport(cl, "run_simulations")


  rep_estimates <- foreach(iter = seq_along(sim), .packages = c("pbapply", "nimble", "MCMCglmm", "coda", "MCMCpack", "boot", "MASS", "parallel","ggmcmc")) %dopar% {
    tryCatch({ run_simulations(sim[[iter]], type = "constant", model_selection = TRUE) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

  }

  #save the results
  save(rep_estimates, file=paste0("constant/estimateDataConstant1",iCount,".RData"))
  stopCluster(cl)
}

