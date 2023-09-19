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
                           inits,
                           returnSamples,
                           nBootReps,
                           parallel = FALSE,
                           silent,
                           IDMmodel){
  message(paste("dropping data fold", i))
  if(parallel) {
    dirName <- file.path(tempdir(), 'nimble_generatedCode', paste0("worker_", i))
  } else dirName = NULL
  model <- conf$model
  leaveOutNames <- model$expandNodeNames(foldFunction(i))
  currentDataNames <- model$getNodeNames(dataOnly = TRUE)
  currentDataNames <- currentDataNames[!(currentDataNames %in% leaveOutNames)]
  saveData <- values(model, leaveOutNames)
  if(!silent) newModel <- model$newModel(check = FALSE, replicate = TRUE)
  else newModel <- suppressMessages(model$newModel(check = FALSE, replicate = TRUE))
  newModel$resetData()
  values(newModel, leaveOutNames) <- NA
  newModel$setData(model$getVarNames(nodes = currentDataNames))
  message(paste("Compiling the model for fold", i))
  if(!silent) Cmodel <-  compileNimble(newModel, dirName = dirName)
  else Cmodel <- suppressMessages(compileNimble(newModel, dirName = dirName))
  predLoss <- FALSE
  if(is.character(lossFunction) && lossFunction == 'predictive'){
    paramNames <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
    dependencies <- model$getDependencies(paramNames)
    paramNamesExpanded <- model$expandNodeNames(paramNames, returnScalarComponents = TRUE)
    #take_out_infinite <-
    #dependencies <- dependencies[model$isData(dependencies)]
    #missingDataNames <- leaveOutNames
    if(!IDMmodel){
      missingDataNames <- leaveOutNames
    }else{
      missingDataNamesX <- leaveOutNames[leaveOutNames %in% model$expandNodeNames(nodes = "X") ]
      missingDataNamesY <- leaveOutNames[leaveOutNames %in% model$expandNodeNames(nodes = "Y") ]
      missingDataNames <- leaveOutNames
      }



    lossFunction <- function(posteriorSamples, observedData){
      start_time <- Sys.time()
      values(model, paramNamesExpanded) <- posteriorSamples
      model$calculate(dependencies)
      if(IDMmodel){
        ret <- model$calculate(missingDataNames)
        ret1 <- model$calculate(missingDataNamesY)
        ret2 <- ret/length(missingDataNames)
        ret3 <- ret/length(missingDataNamesY)
      }else{
        ret <- model$calculate(missingDataNames)
        ret1 <- NA
        ret2 <- ret/length(missingDataNames)
        ret3 <- NA
      }

      end_time <- Sys.time()
      print(time_taken <- end_time - start_time)
      retVal <- c(ret, ret1, ret2, ret3)
      return(retVal)
    }
    leaveOutNames <- paramNames
    predLoss <- TRUE
  }
  message(paste("MCMCconfiguration for fold", i))
  if(!silent) modelMCMCConf <- configureMCMC(Cmodel, nodes = leaveOutNames, monitors = leaveOutNames, print = TRUE, useConjugacy = TRUE)
  #if(!silent) modelMCMCConf <- configureMCMC(Cmodel, monitors = leaveOutNames, print = TRUE, useConjugacy = TRUE)
   else modelMCMCConf <- suppressMessages(configureMCMC(newModel, nodes = leaveOutNames, monitors = leaveOutNames, print = silent, useConjugacy = FALSE))
  if(!predLoss) {
    for(i in seq_along(modelMCMCConf$samplerConfs)) {
      sConf <- modelMCMCConf$samplerConfs[[i]]
      conf$addSampler(target=sConf$target, type=sConf$samplerFunction, control=sConf$control, silent=TRUE)
    }
  }

  # modelMCMCConf$removeSamplers(c("lamLatent",
  #                           "betaSpecies", "betaLatitude", "alphaSpecies", "alphaSites", "betaVisits"), print = FALSE)
  # #mcmcconf$addSampler(target = "delta", type = "RW_block", print = TRUE)
  # #mcmcconf$addSampler(target = "sig", type = "RW_block", print = TRUE)
  # #mcmcconf$addSampler(target = "phi", type = "RW_block", print = TRUE)
  # modelMCMCConf$addSampler(target = "lamLatent", type = "RW_block", print = TRUE)
  # modelMCMCConf$addSampler(target = "betaSpecies", type = "RW_block", print = TRUE)
  # modelMCMCConf$addSampler(target = "betaLatitude", type = "RW_block", print = TRUE)
  # modelMCMCConf$addSampler(target = "alphaSpecies", type = "RW_block", print = TRUE)
  # modelMCMCConf$addSampler(target = "alphaSites", type = "RW_block", print = TRUE)

  MCMCwarnUnsampledStochasticNodes_current <- nimbleOptions('MCMCwarnUnsampledStochasticNodes')
  nimbleOptions(MCMCwarnUnsampledStochasticNodes = FALSE)
  nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)

  if(!silent){
    modelMCMC <- buildMCMC(modelMCMCConf, useConjugacy = TRUE)
    C.modelMCMC <- compileNimble(modelMCMC,
                                 project = Cmodel,
                                 resetFunctions = TRUE,
                                 dirName = dirName)
    C.runmodelMCMC <- runMCMC(C.modelMCMC,
                      niter = niter, #300000
                      nburnin = nburnin,#100000
                      inits = inits,
                      nchains = 1,
                     # thin = 1,
                      summary=TRUE,
                      samples=TRUE,
                      #setSeed = TRUE,
                      #samplesAsCodaMCMC=TRUE,
                      WAIC = FALSE)
  }else{
    modelMCMC <- suppressMessages(buildMCMC(modelMCMCConf))
    C.modelMCMC <- suppressMessages(compileNimble(modelMCMC,
                                                  project = newModel, dirName = dirName))
    silentNull <- suppressMessages(C.modelMCMC$run(niter, progressBar = FALSE))
  }
  leaveOutNamesExpanded <- newModel$expandNodeNames(leaveOutNames, returnScalarComponents = TRUE, sort = TRUE)
  MCMCout <- as.matrix(C.runmodelMCMC$samples)[,leaveOutNamesExpanded, drop = FALSE]
  sampNum <- dim(MCMCout)[1]
  #startIndex <- nburnin+1
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
  #crossVal<- mean(apply(MCMCout[1:sampNum,, drop = FALSE], 1, lossFunction, saveData))
  crossVal<- rowMeans(apply(MCMCout[1:sampNum,, drop = FALSE], 1, lossFunction, saveData))
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

RMSElossFunction <- function(simulatedDataValues, actualDataValues){
  RMSE <- sqrt(mean((simulatedDataValues - actualDataValues)^2))
  return(RMSE)
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
                               lossFunction = 'MSE',
                               MCMCcontrol = list(),
                               returnSamples = FALSE,
                               nCores = 1,
                               nBootReps = 200,
                               silent = FALSE,
                               IDMmodel){

  model <- MCMCconfiguration$model
  niter <- MCMCcontrol[['niter']]
  inits <- MCMCcontrol[['inits']]
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

  if(is.character(lossFunction) && lossFunction == 'RMSE'){
    lossFunction <- RMSElossFunction
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
                                      inits,
                                      returnSamples,
                                      nBootReps,
                                      TRUE,
                                      silent,
                                      mc.cores = nCores,
                                      IDMmodel)
  } else{
    crossValOut <- lapply(1:k, mycalcCrossVal,
                          MCMCconfiguration,
                          foldFunction,
                          lossFunction,
                          niter,
                          nburnin,
                          inits,
                          returnSamples,
                          nBootReps,
                          FALSE,
                          silent,
                          IDMmodel)
  }
  meanCV <- sapply(crossValOut, function(x) x$crossValAverage)
  CVvalue <- rowMeans(meanCV,
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

