
#' Fitting INLA within NIMBLE for Data Generating process
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param data, code, family,n.iterations, n.chains, n.burnin
#'
#' @return MCMC output
#' @export
INLAWiNimDataGeneratingTargetDivide <- function(data,
                                    covariate,
                                    code,
                                    family,
                                    modelData,
                                    modelConstants,
                                    modelInits,
                                    nimbleINLA,
                                    inlaMCMC = c("inla", "mcmc", "inlamcmc"),
                                    inlaMCsampler = "RW_INLA_block",
                                    samplerControl = list(),
                                    parametersToMonitor = list(mcmc = c("mcmc"),
                                                               mcmc2 = c("mcmc"),
                                                               inla = c("inla"),
                                                               additionalPars = NULL),
                                    mcmcConfiguration = list(n.chains = 1,
                                                             n.iterations = 10,
                                                             n.burnin = 0,
                                                             n.thin = 1,
                                                             setSeed = TRUE,
                                                             samples=TRUE,
                                                             samplesAsCodaMCMC = TRUE,
                                                             summary = TRUE,
                                                             WAIC = FALSE)){
  #extract necessary variables
  x <- covariate # must be a matrix
  y <- data # must be a vector
  family <- family
  fixedVals <- parametersToMonitor$inla
  target <- parametersToMonitor$mcmc
  targetMCMC <- parametersToMonitor$mcmcINLA

  #set up initial value

  initsList <- modelInits()

  #Create the model in nimble
  mwtc <- nimble::nimbleModel(code,
                              data = modelData,
                              constants = modelConstants,
                              inits = initsList)

  # Create the model in C
  Cmwtc <- nimble::compileNimble(mwtc,
                                 showCompilerOutput = FALSE) #Have issues compiling

  #create list to return
  retList <- list()
  # Fit INLAMCMC

  if(inlaMCMC %in% "inlamcmc"){
    mcmcconf <- nimble::configureMCMC(Cmwtc,
                                      nodes = c(target))

    # setting sampler controls
    samplerControl$fit.inla = nimbleINLA
    samplerControl$x = x
    samplerControl$y = y
    samplerControl$fixedVals = fixedVals
    samplerControl$fam = family
    #samplerControl$targetMCMC = targetMCMC

    # mcmc configuration
    mcmcconf$addSampler(target = targetMCMC,
                        type = inlaMCsampler,
                        control = samplerControl)
    #mcmcconf$addSampler(target = target[1],
     #                   type = "RW_block")

    mcmcconf$printSamplers(executionOrder = TRUE)
    mcmcconf$addMonitors(c(target, targetMCMC))

    if(!is.null(parametersToMonitor$additionalPars)){
      mcmcconf$addMonitors(parametersToMonitor$additionalPars)
    }

    #build model
    Rmcmc <- nimble::buildMCMC(mcmcconf)

    # Compile
    cmcmc <- nimble::compileNimble(Rmcmc,
                                   project = Cmwtc)
  startTime <- Sys.time()
    mcmc.out <- nimble::runMCMC(cmcmc,
                                niter = mcmcConfiguration[["n.iterations"]],
                                nchains = mcmcConfiguration[["n.chains"]],
                                nburnin = mcmcConfiguration[["n.burnin"]],
                                #inits = initsList,
                                thin = mcmcConfiguration[["n.thin"]],
                                setSeed = mcmcConfiguration[["setSeed"]],
                                samples = mcmcConfiguration[["samples"]],
                                samplesAsCodaMCMC = mcmcConfiguration[["samplesAsCodaMCMC"]],
                                summary = mcmcConfiguration[["summary"]],
                                WAIC = mcmcConfiguration[["WAIC"]])
    endTime <- Sys.time()
    timeTaken <- difftime(endTime, startTime, units = "secs")
    #as.numeric(endTime - startTime)

    ret <- list(mcmc.out = mcmc.out,
                timeTaken = timeTaken)

    #save inlamcmc results
    retList$inlamcmc <- ret
  }

  if(inlaMCMC %in% "mcmc"){
    mcmcconf <- nimble::configureMCMC(Cmwtc,
                                      monitors = c(target, fixedVals, targetMCMC))

    mcmcconf$printSamplers()
    if(!is.null(parametersToMonitor$additionalPars)){
      mcmcconf$addMonitors(parametersToMonitor$additionalPars)
    }
    #build model
    Rmcmc <- nimble::buildMCMC(mcmcconf)

    # Compile
    cmcmc <- nimble::compileNimble(Rmcmc,
                                   project = Cmwtc,
                                   resetFunctions = TRUE)
    startTime <- Sys.time()
    mcmc.out <- nimble::runMCMC(cmcmc,
                                niter = mcmcConfiguration[["n.iterations"]],
                                nchains = mcmcConfiguration[["n.chains"]],
                                nburnin = mcmcConfiguration[["n.burnin"]],
                                #inits = initsList,
                                thin = mcmcConfiguration[["n.thin"]],
                                setSeed = mcmcConfiguration[["setSeed"]],
                                samples = mcmcConfiguration[["samples"]],
                                samplesAsCodaMCMC = mcmcConfiguration[["samplesAsCodaMCMC"]],
                                summary = mcmcConfiguration[["summary"]],
                                WAIC = mcmcConfiguration[["WAIC"]])
    endTime <- Sys.time()
    timeTaken <- difftime(endTime, startTime, units = "secs")

    ret <- list(mcmc.out = mcmc.out,
                timeTaken = timeTaken)

    #save results for MCMC
    retList$mcmc <-  ret
  }
  return(retList)
}


#create matrix

