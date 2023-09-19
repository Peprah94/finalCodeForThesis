# source("crossValidation/cross_validation.R")
# load("CaseStudy/data_idm.RData")
# methodVars <- c("IG", "IG" ,"Spe", "IDM", "IDM")
# modelVars <- c("shared", "covariate","shared","shared","covariate")
# estimatesCaseStudy <- list()
# for(iCount in 1:5){
#   library(doParallel)
#   library(foreach)
#   library(parallel)
#
#   #source("crossValidation/cross_validation.R")
#
#   #cl <- parallel::makeCluster(4)
#   #doParallel::registerDoParallel(cl)
#   #setDefaultCluster(cl)
#
#   #clusterExport(cl, c("run_nimble_model",
#   #                    "mysum",
#   #                    "nimble_sum",
#   #                   "calcCrossValSD",
#   #                    "mycalcCrossVal",
#   #                   "myrunCrossValidate",
#   #                   "MSElossFunction"))
#
#   #estimatesCaseStudy <- list()
#
#   #estimatesCaseStudy <- foreach(iter = seq_along(modelVars), .packages = c("pbapply", "nimble", "MCMCglmm", "coda", "MCMCpack", "boot", "MASS", "parallel","ggmcmc")) %dopar% {
#   estimatesCaseStudy[[iCount]] <- tryCatch({ run_nimble_model(simulations_all[[3]], model = methodVars[iCount],covariance_prior = "LV", method = modelVars[iCount]) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#   #}
#
#
#   #stopCluster(cl)
#
#   #run_nimble_model(simulations_all[[1]], model = methodVars[3],covariance_prior = "LV", method = modelVars[3])
# }

source("crossValidation/cross_validation.R")
load("CaseStudy/data_idm.RData")
methodVars <- c("IG", "IG","Spe", "IDM", "IDM")
modelVars <- c("shared", "covariate", "shared","shared","covariate")
#estimatesCaseStudy <- list()
#for(iCount in 1:5){
library(doParallel)
library(foreach)
library(parallel)


cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)
setDefaultCluster(cl)

clusterExport(cl, c("run_nimble_model",
                    "mysum",
                    # "formatMatrix",
                    "nimble_sum",
                    #"nimbleFormatMatrix",
                    "calcCrossValSD",
                    "mycalcCrossVal",
                    "myrunCrossValidate",
                    "MSElossFunction"))

estimatesCaseStudy <- foreach(iter = seq_along(modelVars), .packages = c("pbapply", "nimble", "MCMCglmm", "coda", "parallel", "foreach", "doParallel")) %dopar% {
  tryCatch({ run_nimble_model(simulations_all[[3]], model = methodVars[iter],covariance_prior = "LV", method = modelVars[iter]) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

  save(estimatesCaseStudy, file=paste0("crossValidation/solitarybees/estimateCrossValidate.RData"))


#}
