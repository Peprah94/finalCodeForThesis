# This script runs the "main" model with all the simulated data
# Takes as input the simulated data, and run_simulations function


library(doParallel)

for(iCount in 1:16){
  message(paste("Running code iteration:", iCount))
load(paste0("simData/simulatedData",iCount,".RData"))
source("all_sim_new.R")

sim <- sim[1:600]

cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
setDefaultCluster(cl)

clusterExport(cl, "run_simulations")


rep_estimates <- foreach(iter = seq_along(sim), .packages = c("pbapply", "nimble", "MCMCglmm", "coda", "MCMCpack", "boot", "MASS", "parallel","ggmcmc")) %dopar% {
  tryCatch({ run_simulations(sim[[iter]], type = "main", model_selection = TRUE) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


save(rep_estimates, file=paste0("main/estimateDataMain1",iCount,".RData"))
stopCluster(cl)
}






