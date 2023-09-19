source("CaseStudy/nimbleCaseStudy.R")

library(doParallel)
library(foreach)
library(parallel)
#methodVars <- c("IDM", "Spe", "IDM")
#modelVars <- c("shared", "shared", "covariate")
methodVars <- c("IG", "IG")
modelVars <- c("shared", "covariate")
load("CaseStudy/data_idm.RData")
source("CaseStudy/nimbleCaseStudy.R")

  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  setDefaultCluster(cl)

  clusterExport(cl, c("run_nimble_model",
                      "mysum",
                      "formatMatrix",
                      "nimble_sum",
                      "nimbleFormatMatrix"))

  estimatesCaseStudy <- foreach(iter = seq_along(modelVars), .packages = c("pbapply", "nimble", "MCMCglmm", "coda", "MCMCpack", "boot", "MASS", "parallel","ggmcmc")) %dopar% {
    tryCatch({ run_nimble_model(simulations_all[[2]], model = methodVars[iter],covariance_prior = "LV", method = modelVars[iter]) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }

  save(estimatesCaseStudy, file = "CaseStudy/hoverflies/estimateCaseStudyIG.RData")
  stopCluster(cl)

#}
