calcCrossValSD <- function(MCMCout, sampNum, nBootReps, saveData, lossFunction, predLoss){
  l <- ceiling(min(1000, (sampNum)/20)) ##length of each
  q <- sampNum - l + 1
  h <- ceiling(sampNum/l)
  blockcvValues <- numeric(nBootReps) ## block estimates of cv
  for(r in 1:nBootReps){
    randNum <- runif(h, 0, 1)
    randIndexStart <- ceiling(randNum*q)
    indexes <- cbind(randIndexStart, randIndexStart + l - 1)
    chunks <- c(apply(indexes, 1, function(x){
      seq(x[1], x[2], 1)}
    ))
    blockcvValues[r] <- mean(apply(MCMCout[chunks,, drop = FALSE], 1, lossFunction, saveData))
    if(predLoss){
      blockcvValues[r] <- log(blockcvValues[r])
    }
  }
  return(sd(blockcvValues))
}


mycalcCrossVal <- function(i,
                         conf,
                         foldFunction,
                         lossFunction,
                         niter,
                         nburnin,
                         returnSamples,
                         nBootReps,
                         parallel = FALSE,
                         silent){
  message(paste("dropping data fold", i))
  if(parallel) {
    dirName <- file.path(tempdir(), 'nimble_generatedCode', paste0("worker_", i))
  } else dirName = NULL
  model <- conf$model
  leaveOutNames <- model$expandNodeNames(foldFunction(i))
  currentDataNames <- model$getNodeNames(dataOnly = TRUE)
  currentDataNames <- currentDataNames[!(currentDataNames %in% leaveOutNames)]
  saveData <- nimble::values(model, leaveOutNames)
  if(!silent) newModel <- model$newModel(check = FALSE, replicate = TRUE)
  else newModel <- suppressMessages(model$newModel(check = FALSE, replicate = TRUE))
  newModel$resetData()
  values(newModel, leaveOutNames) <- NA
  newModel$setData(model$getVarNames(nodes = currentDataNames))
  if(!silent) compileNimble(newModel, dirName = dirName)
  else Cmodel <- suppressMessages(compileNimble(newModel, dirName = dirName))
  predLoss <- FALSE
  if(is.character(lossFunction) && lossFunction == 'predictive'){
    paramNames <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
    dependencies <- model$getDependencies(paramNames)
    missingDataNames <- leaveOutNames
    lossFunction <- function(posteriorSamples, observedData){
      values(model, paramNames) <- posteriorSamples
      model$calculate(dependencies)
      return(model$calculate(missingDataNames))
    }
    leaveOutNames <- paramNames
    predLoss <- TRUE
  }
  if(!silent) modelMCMCConf <- configureMCMC(newModel, nodes = leaveOutNames, monitors = leaveOutNames, print = silent)
  else modelMCMCConf <- suppressMessages(configureMCMC(newModel, nodes = leaveOutNames, monitors = leaveOutNames, print = silent))
  if(!predLoss) {
    for(i in seq_along(modelMCMCConf$samplerConfs)) {
      sConf <- modelMCMCConf$samplerConfs[[i]]
      conf$addSampler(target=sConf$target, type=sConf$samplerFunction, control=sConf$control, silent=TRUE)
    }
  }
  if(!silent){
    modelMCMC <- buildMCMC(modelMCMCConf)
    C.modelMCMC <- compileNimble(modelMCMC,
                                 project = newModel,
                                 dirName = dirName)
    C.modelMCMC$run(niter)
  }
  else{
    modelMCMC <- suppressMessages(buildMCMC(modelMCMCConf))
    C.modelMCMC <- suppressMessages(compileNimble(modelMCMC,
                                                  project = newModel, dirName = dirName))
    silentNull <- suppressMessages(C.modelMCMC$run(niter, progressBar = FALSE))
  }
  leaveOutNamesExpanded <- newModel$expandNodeNames(leaveOutNames, returnScalarComponents = TRUE)
  MCMCout <- as.matrix(C.modelMCMC$mvSamples)[,leaveOutNamesExpanded, drop = FALSE]
  sampNum <- dim(MCMCout)[1]
  startIndex <- nburnin+1
  if(predLoss){
    values(model, missingDataNames) <- saveData
  }
  if(!is.na(nBootReps)){
    crossValAveSD <- calcCrossValSD(MCMCout=MCMCout, sampNum=sampNum,
                                    nBootReps=nBootReps,
                                    saveData=saveData,
                                    lossFunction = lossFunction,
                                    predLoss = predLoss)
  }
  else{
    crossValAveSD <- NA
  }
  crossVal <- mean(apply(MCMCout[startIndex:sampNum,, drop = FALSE], 1, lossFunction, saveData))
  if(predLoss){
    crossVal <- crossVal
  }
  return(list(crossValAverage = crossVal,
              crossValAveSD = crossValAveSD,
              samples= if(returnSamples) MCMCout else NA))
}

MSElossFunction <- function(simulatedDataValues, actualDataValues){
  MSE <- mean((simulatedDataValues - actualDataValues)^2)
  return(MSE)
}

generateRandomFoldFunction <- function(model, k){
  dataNodes <- model$expandNodeNames(model$getNodeNames(dataOnly = TRUE))
  nDataNodes <- length(dataNodes)
  if(k > nDataNodes){
    stop('Cannot have more folds than there are data nodes in the model.')
  }
  approxNumPerFold <- floor(nDataNodes/k)
  leaveOutNames <- list()
  reducedDataNodeNames <- dataNodes
  for(i in 1:(k-1)){
    leaveOutNames[[i]] <- sample(reducedDataNodeNames, approxNumPerFold)
    reducedDataNodeNames <- reducedDataNodeNames[-which(reducedDataNodeNames %in% leaveOutNames[[i]])]
    approxNumPerFold <- floor(length(reducedDataNodeNames)/(k-i))
  }
  leaveOutNames[[k]] <- reducedDataNodeNames
  randomFoldFunction <- function(i){
    return(leaveOutNames[[i]])
  }
}

 
myrunCrossValidate <- function(MCMCconfiguration,
                             k,
                             foldFunction = 'random',
                             lossFunction = 'predictive',
                             MCMCcontrol = list(), 
                             returnSamples = FALSE,
                             nCores = 1,
                             nBootReps = 200,
                             silent = FALSE){
  
  model <- MCMCconfiguration$model
  niter <- MCMCcontrol[['niter']] 
  if(k < 2){
    stop("k must be at least 2.")
  }
  if(is.null(niter)){
    niter <- 10000
    warning("Defaulting to 10,000 MCMC iterations for each MCMC run.")
  }
  nburnin <-  MCMCcontrol[['nburnin']]  
  if(is.null(nburnin)){
    nburnin <- floor(0.1*niter)
  }else{if(nburnin >= niter | nburnin < 0){
    stop("nburnin needs to be between 0 and niter.")
  }}
  if(is.character(foldFunction) && foldFunction == 'random'){
    foldFunction <- generateRandomFoldFunction(model, k)
  }
  if(is.character(lossFunction) && lossFunction == 'MSE'){
    lossFunction <- MSElossFunction
  }
  allLeaveoutDataNames <- lapply(1:k, function(x){
    thisDataNames <- try(model$expandNodeNames(foldFunction(x)), silent = TRUE)
    if(inherits(thisDataNames, "try-error")){
      stop(paste("foldFunction is returning incorrect node names for fold", x))
    }
    else{
      if(!all(model$expandNodeNames(thisDataNames) %in% model$expandNodeNames(model$getNodeNames(dataOnly = TRUE)))){
        stop(paste("foldFunction is returning names of non-data nodes for fold", x))
      }
    }
  })
  
  if(!requireNamespace('parallel', quietly = TRUE)) {
    warning("To use multiple cores for the cross-validation calculation, the 'parallel' package must be available.  Defaulting to one core.")
    nCores <- 1
  }
  
  if(nCores > 1){
    crossValOut <- parallel::mclapply(1:k, mycalcCrossVal,
                                      MCMCconfiguration,
                                      foldFunction,
                                      lossFunction,
                                      niter,
                                      nburnin,
                                      returnSamples,
                                      nBootReps,
                                      TRUE,
                                      silent,
                                      mc.cores = nCores)
  } else{
    crossValOut <- lapply(1:k, mycalcCrossVal,
                          MCMCconfiguration,
                          foldFunction,
                          lossFunction,
                          niter,
                          nburnin,
                          returnSamples,
                          nBootReps,
                          FALSE,
                          silent)
  }
  meanCV <- sapply(crossValOut, function(x) x$crossValAverage)
  CVvalue <- mean(meanCV[!is.infinite(meanCV)],
                  na.rm=TRUE)
  if(!is.na(nBootReps)){
    CVstandardError <- sqrt(sum(sapply(crossValOut,
                                       function(x) x$crossValAveSD^2),
                                na.rm=TRUE)/k^2)
  }
  else{
    CVstandardError <- NA
  }
  foldCVinfo <- lapply(crossValOut, function(x){return(c(foldCVvalue = x$crossValAverage,
                                                         foldCVstandardError = x$crossValAveSD))})
  out <- list(CVvalue=CVvalue,
              CVstandardError=CVstandardError,
              foldCVinfo = foldCVinfo)
  if(!returnSamples){
    out$samples <- NULL
  }
  else{
    out$samples <- lapply(crossValOut, function(x) x$samples)
  }
  return(out)
}

# myrunCrossValidate(MCMCconfiguration = mcmcconf, 
#                  k = 2,
#                  nCores = 1,
#                  lossFunction = "predictive",
#                  MCMCcontrol = list(niter = 10000, #10000
#                                     nburnin = 9990), #5000
#                  foldFunction = "random",
#                  nBootReps = NA
# )
