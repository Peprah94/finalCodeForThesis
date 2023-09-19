setwd("/Volumes/kwakupa/misclassification/new")
load("/Volumes/kwakupa/misclassification/new/nimble_data.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_constant2.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_intercept2.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_main2.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_only_principal_cov2.RData")
#only_principal_cov <- variable
load("/Volumes/kwakupa/misclassification/new/estimated_data_variable2.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_op_cov_finter2.RData")
load("/Volumes/kwakupa/misclassification/new/estimated_data_ml_estimates2.RData")


#experience <- unique(as.numeric(as.factor(dataout$no_previous_obs)))
experience <- ifelse(dataout$no_previous_obs + 1 > 9, 10,  dataout$no_previous_obs + 1)
nusers <- length(unique(experience))

all_data <- list(variable, 
                 constant, 
                 intercept,
                 main, 
                 op_cov_finter,
                 only_principal_cov,
                 ml_estimates)
#packages
library(ggplot2)
library(dplyr)
theme_set(theme_classic(base_size = 10))


subset_parameters <- function(data, subset_index){
  ret <-  lapply(data, function(x){
    y <- x[[1]]
    if(subset_index == "p"){
      rr <- y[rownames(y)[grepl(subset_index, rownames(y))], 2]
      rr[1]
    }else{
      y[rownames(y)[grepl(subset_index, rownames(y))], 2]
    }
  })
  return(ret)
}

gamma0 <- subset_parameters(all_data, "gamma0")
gamma1 <- subset_parameters(all_data,"gamma1")
#user_effect <- subset_parameters_annex(all_data, 12)

classifyProbs <- lapply(all_data[1:6], function(x){

  nusers
  tr <- lapply(1:nusers, function(i){
    c(paste0("z.omega[1, 1, ", i, "]", collapse = ","),
      paste0("z.omega[2, 2, ", i, "]", collapse = ","),
      paste0("z.omega[3, 3, ", i, "]", collapse = ","),
    paste0("z.omega[4, 4, ", i, "]", collapse = ",")#,
  #paste0("z.omega[5, 5, ", i, "]", collapse = ","
         #)
  )
  })%>%
    unlist()
  
  tt <- sum(!is.na(x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 2]))
  if(tt != 200){
  main <- rep(x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 2][c(1,6,11, 16)],nusers)
 lower <- rep(x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 4][c(1,6,11, 16)],nusers)
 upper <- rep(x$output[rownames(x$output)[grepl("z.omega", rownames(x$output))], 5][c(1,6,11, 16)],nusers)
  }else{
    main <- x$output[tr, 2]
    lower <- x$output[tr, 4]
    upper <- x$output[tr, 5] 
  }
 N <- rep(1:nusers, each = 4)
 #i  <- i +1
 #print(i)
 ret <- cbind(N, main, lower, upper)
 return(ret)
})%>%
  do.call("rbind", .)%>%
  data.frame()%>%
  dplyr::mutate(Method = c(rep("covariate", 4*nusers),
                    rep("constant", 4*nusers),
                    rep("intercept", 4*nusers),
                    rep("main", 4*nusers),
                    rep("fixed_intercov", 4*nusers),
                    rep("fixed_covariate", 4*nusers)),
         Species = rep(c("i) L. canus",
                         "ii) L. marinus",
                         "iii) L. argentatus",
                         "iv) L. fuscus"
                          ),6*nusers))%>%
  ggplot()+
  #geom_point(aes(x = N, y = value, col = Method))+
  geom_line(aes(x = as.factor(N), y = main, col = Method, group = Method))+
  geom_ribbon(aes(x = N, ymin= lower, ymax= upper, fill= Method, group = Method), alpha = 0.1)+
  ylab("Probability of correct classification")+
  xlab("")+
  #scale_x_discrete("N", labels = as.character(N), breaks = N)
  facet_wrap(~Species, ncol = 2, scales = "free")

experiencePlots <- ggplot(data = data.frame(experience = experience))+
  geom_bar(aes(x = as.factor(experience)))+
  xlab("experience")

fig <- ggpubr::ggarrange(classifyProbs, 
                         experiencePlots, 
                         nrow =2 , 
                         ncol = 1, 
                         legend = "bottom", heights = c(3,1),
                         labels = c("A", "B"), common.legend = TRUE)

ggsave("classifyProbs.png", 
       plot = fig, 
       dpi = 320, 
       height = 16, 
       width = 14, 
       units = "cm")







#all_data <- list(variable, intercept, constant, main)

accuracy <- subset_parameters(all_data, "accuracy")%>%
  unlist()

precision <- subset_parameters(all_data, "precision")%>%
  unlist()

recall <- subset_parameters(all_data, "recall")%>%
  unlist()

varSelectionProb <- subset_parameters(all_data, "p")%>%
  unlist()

waic <- unlist(lapply(all_data, function(x){
  x$waic$WAIC
}))

method = c("covariate", "constant", "intercept", "main","fixed-intercov", "fixed-covariate", "machine_learning")

predPerformance <- cbind(method, 
                         accuracy, 
                         precision, 
                         recall , 
                         varSelectionProb, 
                         waic)
write.csv(predPerformance,"predPerformanceGullDataset.csv", row.names = FALSE)



beta0 <- subset_parameters(all_data[1:6], "beta0")%>%
  unlist()#%>%
  #matrix(., ncol = 3, byrow = TRUE)

beta1 <- subset_parameters(all_data[1:6], "beta1")%>%
  unlist()

beta11 <- lapply(as.list(1:3), function(i){
  x <- paste0("beta1[", i, ", 1]")
   beta1[names(beta1) == x]
  })%>%
  do.call("rbind", .)%>%
  c()

beta12 <- lapply(as.list(1:3), function(i){
  x <- paste0("beta1[", i, ", 2]")
  beta1[names(beta1) == x]
})%>%
  do.call("rbind", .)%>%
  c()

beta5 <- subset_parameters(all_data[1:6], "beta5")%>%
  unlist()

parameters_estimates <- cbind(beta0, beta11, beta5)%>%
  data.frame()%>%
 dplyr::mutate(method = c(rep("covariate", 3),
                          rep("constant", 3),
                          rep("intercept", 3),
                          rep("main", 3),
                          rep("fixed-intercov", 3),
                          rep("fixed-covariate", 3)),
               species = rep(c("L. canus",
                           "L. marinus",
                           "L. argentatus"), 6))%>%
  reshape2::melt(id.vars = c("method", "species"))%>%
  reshape2::dcast(., species+variable ~ method)%>%
  dplyr::arrange(variable)%>%
  dplyr::mutate(parameter = rep(c("Intercept", 
                                   "Precipitation",
                                  "Experience"), each = 3))
write.csv(parameters_estimates,"ecologicalParametersEstimatesGullData.csv", row.names = FALSE)



###############
#   Classification process paramaters
##############
gamma0 <- subset_parameters(all_data[1:6], "gamma0")%>%
  unlist()


gamma01 <- lapply(as.list(1:4), function(i){
  x <- paste0("gamma0[", i, ", 1]")
  gamma0[names(gamma0) == x]
})%>%
  do.call("rbind", .)%>%
  c()

gamma02 <- lapply(as.list(1:4), function(i){
  x <- paste0("gamma0[", i, ", 2]")
  gamma0[names(gamma0) == x]
})%>%
  do.call("rbind", .)%>%
  c()

gamma03 <- lapply(as.list(1:4), function(i){
  x <- paste0("gamma0[", i, ", 2]")
  gamma0[names(gamma0) == x]
})%>%
  do.call("rbind", .)%>%
  c()

gamma04 <- lapply(as.list(1:4), function(i){
  x <- paste0("gamma0[", i, ", 4]")
  gamma0[names(gamma0) == x]
})%>%
  do.call("rbind", .)%>%
  c()

allgamma0 <- cbind(gamma01, gamma02, gamma03, gamma04)


#gamma1 <- subset_parameters_annex(all_data, 8)
gamma1 <- subset_parameters(all_data[1:6], "gamma1")%>%
  unlist()


gamma11 <- lapply(as.list(1:4), function(i){
  x <- paste0("gamma1[", i, ", 1]")
  gamma1[names(gamma1) == x]
})%>%
  do.call("rbind", .)%>%
  c()

gamma12 <- lapply(as.list(1:4), function(i){
  x <- paste0("gamma1[", i, ", 2]")
  gamma1[names(gamma1) == x]
})%>%
  do.call("rbind", .)%>%
  c()

gamma13 <- lapply(as.list(1:4), function(i){
  x <- paste0("gamma1[", i, ", 2]")
  gamma1[names(gamma1) == x]
})%>%
  do.call("rbind", .)%>%
  c()

gamma14 <- lapply(as.list(1:4), function(i){
  x <- paste0("gamma1[", i, ", 4]")
  gamma1[names(gamma1) == x]
})%>%
  do.call("rbind", .)%>%
  c()

allgamma1 <- cbind(gamma11, gamma12, gamma13, gamma14)

missclassPars <- cbind(allgamma0, allgamma1)%>%
  data.frame()%>%
  dplyr::mutate(method = c(rep("covariate", 4),
                           rep("constant", 4),
                           rep("intercept", 4),
                           rep("main", 4),
                           rep("fixed-intercov", 4),
                           rep("fixed-covariate", 4)),
                species = rep(c("L. canus",
                                "L. marinus",
                                "L. argentatus",
                                "Larus fuscus"), 6))%>%
  reshape2::melt(id.vars = c("method", "species"))%>%
  dplyr::mutate(reported_states = rep(c("L. canus",
                                      "L. marinus",
                                      "L. argentatus",
                                      "Larus fuscus", "L. canus",
                                      "L. marinus",
                                      "L. argentatus",
                                      "Larus fuscus"), each = 24))%>%
  reshape2::dcast(., species + reported_states + variable ~ method)%>%
  dplyr::arrange( reported_states, variable)


write.csv(missclassPars,"classificationParametersEstimatesGullData.csv", row.names = FALSE)






