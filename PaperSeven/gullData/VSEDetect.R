#load data
library(INLA)
library(inlabru)
library(dplyr)
library(sf)
library(terra)
library(raster)
library(dplyr)
library(fmesher)
bru_options_set(inla.mode = "experimental")
#bru_safe_sp(force = TRUE)
bru_options_set(control.compute = list(dic = TRUE,
                                       cpo = TRUE,
                                       waic = TRUE))

load("allDataForModelNew.RData")

#Extract parameters needed
source("gullAnalysis.R")


#ips <- fm_int(mesh, samplers = aa)
# formula = as.formula(paste0("coordinates ~ betaPO0",i,"  + cov1",
#                             i,"  + cov4",
#                             i," + w1",i,"+cov3",i, "+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)+ beta0det",
#                             i,"+timeSpent",i,"+fun1(beta0det",i,", timeSpent",i,") + log(weight)")),
#Likelihoods of INLABru
likPO <- likPA <- likNina <- likNibc <- likDet <- list()
hyperNB <- list()
#Estimates of NB size from fitted model to NINA dataset
k.pm <- c(1.1197, 1.448)
for(i in 1:nspecies){
  hyperNB[[i]] <- list(size = list(initial = k.pm[i],
                                   fixed = TRUE))  
}

for(i in 1:nspecies){
  likPO[[i]] <- inlabru::like("cp",
                              formula = as.formula(paste0("geometry ~ betaPo0",i,"+w1",i, "+cov1", i, "+cov3", i, "+cov4", i, "+ cov2", i)),
                              data = dataBeforeClassification[[i]],
                              #components = cmp,
                              domain = list(geometry = mesh),
                              #ips = fm_int(mesh, samplers = dataBeforeClassification[[i]]),
                              samplers = aa
  )
  # lik2[[i]] <- inlabru::like("cp",
  #                            formula = coordinates ~ beta0thin + cov2 + w2,
  #                            data = spPointsGbif,
  #                            #components = cmp,
  #                            domain = list(coordinates = mesh),
  #                            #ips = ips#,
  #                            samplers = aa
  # )
  # likPA[[i]] <- inlabru::like("binomial",
  #                             formula = as.formula(paste0("detdata",i," ~ betaPA0",i,"+ cov1",
  #                                                         i,"+cov3",i, "+cov4",i,"+w1",i)),
  #                             data = occurenceDf[[i]],
  #                             Ntrials = rep(1,length(detdata[[i]]$Ntrials)),
  #                             control.family= list(link = 'cloglog')#,
  #                             #ips = fm_int(mesh, samplers = detdata[[1]])#,
  #                             #components = cmp,
  #                             #domain = list(geometry = mesh),
  #                             #samplers = aa
  # )
  
  likPA[[i]] <- inlabru::like("binomial",
                              formula = as.formula(paste0("detdata",i," ~ beta0det",i," + timeSpent",i,"+betaPA0",i,"+ cov1",
                                                          i,"+cov3",i, "+cov4",i,"+w1",i)),
                              data = detdata[[i]],
                              Ntrials = detdata[[i]]$Ntrials,
                              control.family= list(link = 'cloglog')#,
                              #ips = fm_int(mesh, samplers = detdata[[1]])#,
                              #components = cmp,
                             # domain = list(geometry = mesh),
                              #samplers = aa
  )
  
  likNina[[i]] <- inlabru::like("nbinomial",
                                formula = as.formula(paste0("individualCount ~ betaNina0",i, "+cov1", i, "+cov3", i, "+cov4", i,"+w1",i)),
                                data = countDf[[i]],
                                #Ntrials = detdata[[i]]$Ntrials,
                                control.family= list(hyper = hyperNB[[i]])#,
                                #ips = fm_int(mesh, samplers = countDf[[i]]),
                                #components = cmp,
                                #domain = list(geometry = mesh),
                                #samplers = aa
  )
  
  likNibc[[i]] <- inlabru::like("poisson",
                                formula = as.formula(paste0("individualCount ~ betaNibc0",i,"  + cov1",
                                                            i,"+cov3",i, "+cov4",i,"+w1",i)),
                                data = nbic[[i]],
                                #Ntrials = detdata[[i]]$Ntrials,
                                 control.family= list(link = 'log')#,
                                #ips = fm_int(mesh, samplers = nbic[[i]]),
                                #components = cmp,
                                #domain = list(geometry = mesh),
                                #samplers = aa
  )
}

# elev$cov11 <- elev$cov12  <- elev$alt
# dist$cov2 <- dist$distanceRaster
# prec$cov31 <- prec$cov32 <- prec$prec01
# timeSpent$timeSpent1 <- timeSpent$timeSpent2 <- timeSpent$layer
# distwb$cov41 <- distwb$cov42 <- distwb$distanceRasterWaterBody

fit4 <- inlabru::bru(cmp,
                     likPO[[1]], likPO[[2]],
                     #unquote(paste0("lik1[[", 1:nspecies, "]]", collapse = ",")),
                     #lik2[[1]],
                     #unquote(paste0("lik3[[", 1:nspecies, "]]", collapse = ",")),
                     likPA[[1]],likPA[[2]],
                     #likDet[[1]],likDet[[2]],
                     #likNina[[1]],likNina[[2]],
                     likNibc[[1]], likNibc[[2]],
                     options = list(control.inla = list(int.strategy = "eb",
                                                        control.vb=list(enable=FALSE),
                                                                    cmin = 0.1
                     ),
                     bru_method = list(
                       #taylor = "pandemic",
                       # search = "all",
                       # factor = (1 + sqrt(5)) / 2,
                       rel_tol = 0.1
                       # max_step = 2,
                       # lin_opt_method = "onestep"
                     ), #change to 0.01
                     bru_max_iter = 10,
                     bru_verbose = 2
                     )
)

save(fit4, file = "VSEDetectFit3.RData")

# e.int <- predict(fit4, pixels(mesh, mask = boundary), ~ exp(beta01 + cov31 + cov11+ w11 + beta0thin + cov2 + beta0det1 + timeSpent1))
# 
# ggplot() +
#   gg(e.int) +
#   gg(boundary, alpha = 0) +
#   #gg(nests, shape = "+") +
#   coord_equal()

