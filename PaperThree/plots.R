# Results file sizes are large. We have attached the link to access our results
# https://www.dropbox.com/sh/9gdf426e2sw5yas/AADqlpM6jxvqRQMmccku3WxTa?dl=0

# We have also managed to attach a ZIP file of the results

# Load data
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)


# set directory to the path the data was downloaded

###################
# Linear Gaussian SSM
####################

load("linearGaussianSSM/simulatedDataEx1.RData")

# Need to format the data from Example One
#write a function to do that
formatEx1 <- function(data){

  result <- lapply(data, function(x){
    BMC <- list(x[[1]], x[[2]], x[[3]])
    BSC <- list(x[[4]], x[[5]], x[[6]])
    putTogether <- list(BMC, BSC)
    others <- x[7:26]
    ret <- c(putTogether,
             others)
    return(ret)
  })

  return(result)
}

# Estimate root mean square errr
averageRMSE <- function(x,y){
  ret <- sqrt(mean((x-y)^2))
  return(ret)
}

## Load data for Aux PF
load("linearGaussianSSM/estimates1.RData")
pf1 <- bootstrapEstimates[1:15]
rm(bootstrapEstimates)

load("linearGaussianSSM/estimates2.RData")
pf2 <- auxiliaryEstimates[1:15]
rm(auxiliaryEstimates)

load("linearGaussianSSM/estimates3.RData")
pf3 <- auxiliaryEstimates[1:15]
rm(auxiliaryEstimates)

# Put all dataset together
allPF <- c(pf1, pf2, pf3)
rm(pf1, pf2, pf3)

## Load Bootstrap results
# load("linearGaussianSSM/bootstrapPF/estimatesBFNew1.RData")
# bootPF1 <- bootstrapEstimates[1:15]
# rm(bootstrapEstimates)
#
# load("linearGaussianSSM/bootstrapPF/estimatesBFNew2.RData")
# bootPF2 <- bootstrapEstimates[1:15]
# rm(bootstrapEstimates)
#
# # Put all dataset together
# bootPF <- c(bootPF1, bootPF2)
# rm(bootPF1, bootPF2)

## Put all results together

# results1 <- c(auxPF,
#               bootPF)

# Flatten list
results <- purrr::flatten(allPF)


## Extract model Parameters

nodeNames = c("ahat", "chat")

#estimate the bias of model parameter
biasModelPars <- lapply(results, function(x){
  red <- x[[2]]$all.chains
  output <-  red[rownames(red)[rownames(red) %in% nodeNames],"Median"]
  return(output)
})%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  dplyr::mutate(t = rep(c(rep(50, 3),
                          rep(c(49, 45, 20, 10, 5), 3)), 45),
                model = rep(
                  c(rep(c("BBSC", "ABSC","BMC",
                          rep("RMC", 5),
                          rep("BUMC", 5),
                          rep("AUMC", 5)), 45))
                ))

## Extract results for auxiliary PF and plot
biasModelPars%>%
  ggplot(., aes(x = as.factor(t), y = ahat, fill = model))+
  geom_boxplot()


# Fill colors for plotting
fill.colors <- c("BMC" = "#FF6600",
                 "ABSC" = "#FF3300",
                 "RMC" = "#00CCFF",
                 "BBSC" = "#E69F00",
                 "AUMC" = "#0033FF",
                 "BUMC" = "#3399FF",
                 "truth" = "black")

# extract results
biasModelPars1 <- biasModelPars%>%
  dplyr::group_by(t,model )%>%
  dplyr::summarise(across(c(ahat,chat), mean))%>%
  reshape2::melt(., id.vars = c("t", "model"))%>%
  # filter(model %in% c("ABMC", "ABSC",
  #                     "ARMC",
  #                     "AUMC", "BUMC", "BBSC"))%>%
  dplyr::mutate(model = replace(model, model == "ABMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "ARMC", "RMC"))


#plot results
biasModelParsAhat <- biasModelPars1%>%
  filter(variable %in% c("ahat"))%>%
  filter( !t %in% c("45","50"))%>%
  ggplot(.,
                             mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 8)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("ahat"),4][1], linetype = "dashed", col = "#FF3300", linewidth = 2)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("ahat"),4][2], linetype = "dotdash", col = "#E69F00", linewidth = 2)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("ahat"),4][3], linetype = "dotted", col = "#FF6600", linewidth = 2)+
  facet_wrap(~variable, ncol = 3, scales = "free_y",
             labeller = as_labeller(c( "ahat" = "a",
                                       #"bhat" = "b",
                                       "chat" = "c")))+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors)+
  xlab("")+
  ylab("Bias")

biasModelParsChat <- biasModelPars1%>%
  filter(variable %in% c("chat"))%>%
  filter( !t %in% c("45","50"))%>%
  ggplot(.,
         mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 8)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("chat"),4][1], linetype = "dashed", col = "#FF3300", linewidth = 2)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("chat"),4][2], linetype = "dotdash", col = "#E69F00", linewidth = 2)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("chat"),4][3], linetype = "dotted", col = "#FF6600", linewidth = 2)+
  facet_wrap(~variable, ncol = 3, scales = "free_y",
             labeller = as_labeller(c( "ahat" = "a",
                                       #"bhat" = "b",
                                       "chat" = "c")))+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors)+
  xlab(" ")+
  ylab(" ")


### Extract RMSE results and plot
# Estimating MCSE
mcseEstimates <- lapply(results, function(x){mcmcse::mcse.mat(as.matrix(rbind(x[[1]]$chain1,
                                                                              x[[1]]$chain2))[, c("a", "c")],
                                                              method = "bm",
                                                              #size = bacthSize,
                                                              g = NULL)%>%
    as.data.frame()%>%
    dplyr::select(se)%>%
    t()})%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
dplyr::mutate(t = rep(c(rep(50, 3),
                                                    rep(c(49, 45, 20, 10, 5), 3)), 45),
                                          model = rep(
                                            c(rep(c("BBSC", "ABSC","BMC",
                                                    rep("RMC", 5),
                                                    rep("BUMC", 5),
                                                    rep("AUMC", 5)), 45))
                                          ))%>%
  reshape2::melt(., id.vars = c("t", "model"))%>%
  dplyr::group_by(model, t, variable)%>%
  #dplyr::group_by(t)%>%
  dplyr::summarise(value = mean(value))

#select for auxPF
mcseEstimatesPars1 <- mcseEstimates%>%
  # filter(model %in% c("ABMC", "ABSC",
  #                     "ARMC",
  #                     "AUMC", "BUMC", "BBSC"))%>%
  dplyr::mutate(model = replace(model, model == "ABMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "ARMC", "RMC"))

mcseModelParsA <- mcseEstimatesPars1%>%
  filter( !t %in% c("45","50"))%>%
  filter(variable %in% c("a"))%>%
  ggplot(data = .,
                             mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 8)+
  geom_hline(yintercept = c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("a"),4]$value)[1], linetype = "dashed", col = "#FF3300",linewidth = 2)+
  geom_hline(yintercept = c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("a"),4]$value)[2], linetype = "dotdash", col = "#E69F00",linewidth = 2)+
  geom_hline(yintercept =c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("a"),4]$value)[3], linetype = "dotted", col = "#FF6600",linewidth = 2)+
  facet_wrap(~variable, ncol = 3, scales = "free_y")+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors)+
  xlab(" ")+
  ylab("MCSE")

# Extract bootstrap PF
mcseModelParsB <- mcseEstimatesPars1%>%
  filter( !t %in% c("45","50"))%>%
  filter(variable %in% c("c"))%>%
  ggplot(data = .,
         mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 8)+
  geom_hline(yintercept = c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("c"),4]$value)[1], linetype = "dashed", col = "#FF3300",linewidth = 2)+
  geom_hline(yintercept = c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("c"),4]$value)[2], linetype = "dotdash", col = "#E69F00",linewidth = 2)+
  geom_hline(yintercept =c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("c"),4]$value)[3], linetype = "dotted", col = "#FF6600",linewidth = 2)+
  facet_wrap(~variable, ncol = 3, scales = "free_y")+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors)+
  xlab(" ")+
  ylab(" ")



## Results for latent state distribution
# auxiliary latent states
auxLatentEst <- lapply(seq_along(allPF), function(i){

  x <- allPF[[i]]
  y <- simData[[i]]$x

  BBSC <- averageRMSE(x[[1]][[2]]$all.chains[!rownames(x[[1]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)
  ABSC <- averageRMSE(x[[2]][[2]]$all.chains[!rownames(x[[2]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)
  BMC <- averageRMSE(x[[3]][[2]]$all.chains[!rownames(x[[3]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)

  # BUMC
  BUMC_49 <- averageRMSE(c(x[[9]][[2]]$all.chains[!rownames(x[[9]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  BUMC_45 <- averageRMSE(c(x[[10]][[2]]$all.chains[!rownames(x[[10]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  BUMC_20 <- averageRMSE(c(x[[11]][[2]]$all.chains[!rownames(x[[11]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  BUMC_10 <- averageRMSE(c(x[[12]][[2]]$all.chains[!rownames(x[[12]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  BUMC_5 <- averageRMSE(c(x[[13]][[2]]$all.chains[!rownames(x[[13]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  # AUMC
  AUMC_49 <- averageRMSE(c(x[[14]][[2]]$all.chains[!rownames(x[[14]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  AUMC_45 <- averageRMSE(c(x[[15]][[2]]$all.chains[!rownames(x[[15]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  AUMC_20 <- averageRMSE(c(x[[16]][[2]]$all.chains[!rownames(x[[16]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  AUMC_10 <- averageRMSE(c(x[[17]][[2]]$all.chains[!rownames(x[[17]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  AUMC_5 <- averageRMSE(c(x[[18]][[2]]$all.chains[!rownames(x[[18]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  ret <- data.frame(BBSC, ABSC, BMC,
                    BUMC_49, BUMC_45,
                    BUMC_20, BUMC_10,
                    BUMC_5,AUMC_49, AUMC_45,
                    AUMC_20, AUMC_10,
                    AUMC_5
  )
  return(ret)
})%>%
  do.call("rbind", .)%>%
  reshape2::melt()%>%
  dplyr::group_by(variable)%>%
  dplyr::summarise(mean = mean(value),
                   sd = sd(value))%>%
  dplyr::mutate(t = c(rep(50, 3),
                          rep(c(49, 45, 20, 10, 5), 2)),
                model = c("BBSC", "ABSC","BMC",
                          rep("BUMC", 5),
                          rep("AUMC", 5))
                )%>%
  dplyr::mutate(model = replace(model, model == "ABMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "ARMC", "RMC"))

#Plot error plot
errorModelParsAuxPF <- auxLatentEst%>%
  filter( !t %in% c("45","50"))%>%
 # filter(variable %in% c("c"))%>%
  ggplot(data = ., mapping = aes(x = as.factor(t), y = mean, col = model))+
  geom_point(position = position_dodge(width = 0.7), size = 8)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),size=1)+
  geom_hline(yintercept = c(auxLatentEst[auxLatentEst$t==50,2]$mean)[1], linetype = "dashed", col = "#FF3300",linewidth = 2)+
  geom_hline(yintercept = c(auxLatentEst[auxLatentEst$t==50,2]$mean)[2], linetype = "dotdash", col = "#E69F00",linewidth = 2)+
  geom_hline(yintercept =c(auxLatentEst[auxLatentEst$t==50,2]$mean)[3], linetype = "dotted", col = "#FF6600",linewidth = 2)+
   theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors)+
  xlab( " ")+
  ylab(expression(paste("Mean ", "\u00B1", " SE")))


#Save results
biasPlots <- ggarrange(biasModelParsAhat,
                      biasModelParsChat,
                      mcseModelParsA,
                      mcseModelParsB,
                      errorModelParsAuxPF,
                      nrow = 3,
                      ncol = 2,
                      widths = c(1,1,2),
                     # labels = c("A)", "B)", "C"),
                      font.label = list(size = 25),
                      common.legend = TRUE)

auxPFpLots <- annotate_figure(biasPlots,
                bottom = text_grob(expression(t[r]),
                                   size = 25))


## Put plots together

ggsave(filename = "Figures/auxPFplots.png",
       plot = auxPFpLots,
       width = 22,
       height = 20,
       dpi = 150)

### Convergence plots

## Auxiliary PF
#MCMC
auxPFSubset <- allPF[[1]][c(1,2, 3, 9, 13, 14, 18)]

aPlots <- lapply(auxPFSubset, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% c("a", "c", "x[1]", "x[2]"))%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")+
    theme(legend.position = "bottom")
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots[1:3],
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D", "E"))

aPlots1 <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90),
                          bottom = text_grob("Iterations"))

ggsave(filename = "Figures/baselineModels.png",
       plot = aPlots1,
       width = 6,
       height = 7,
       units = "in")

fig <- ggpubr::ggarrange(plotlist = aPlots[4:7],
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D", "E"))
aPlots2 <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90),
                          bottom = text_grob("Iterations"))

ggsave(filename = "Figures/updatedModels.png",
       plot = aPlots2,
       width = 8,
       height = 8,
       units = "in")



rm(list=ls())




######################
# Simulation study 2: Dynamic occupancy model
###################

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)

# Fill colors for plotting
fill.colors <- c("BMC" = "#FF6600",
                 "ABSC" = "#FF3300",
                 "RMC" = "#00CCFF",
                 "BBSC" = "#E69F00",
                 "AUMC" = "#0033FF",
                 "BUMC" = "#3399FF",
                 "truth" = "black")

#Load simulated data
load("dynamicOccupanyModel/Bootstrap/simDataDynamicOccupancy.RData")

#load bootstrap results
load("dynamicOccupanyModel/Bootstrap/example4BaselineMCMC1.RData")
BBMC <- baselineMCMC

load("dynamicOccupanyModel/Bootstrap/example4ReducedBootstrapTrue.RData")
BRMC <- example2ReducedModelTrue

load("dynamicOccupanyModel/Bootstrap/example4UpdatedBootstrapTrue1.RData")
BUMC <- example2UpdatedModelTrue

load("dynamicOccupanyModel/auxiliaryPF/example2UpdatedAuxiliaryTrue1.RData")
AUMC <- example2UpdatedModelTrue


#remove individual data results
rm(baselineMCMC,
   example2ReducedModelTrue,
   example2UpdatedModelTrue)

allModels <- list(BBMC,BRMC, BUMC, AUMC)


## Extract realised occupancy and plot
ret <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Mean"]
})

# Realised occupancy
occSites <- data.frame(BMC = ret[[1]],
                       RMC = c(ret[[2]], rep(NA, 10)),
                       BUMC = c(ret[[3]]),
                       AUMC = c(ret[[4]]),
                       truth = simData$occSites,
                       year = 1:55)%>%
  reshape2::melt(id.vars = c("year"))%>%
  ggplot(data = ., mapping = aes(x = year,
                                 y = value, col = variable, group = variable))+
  geom_point(position = position_dodge(width = 1),size=5)+
  geom_line(position = position_dodge(width =1),linewidth = 3)+
  theme_classic()+
  #ylim(c(0.75, 1.6))+
  geom_vline(xintercept = 45, linetype = "dashed", col = "red")+
  scale_color_manual(name = "Model", values = fill.colors)+
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50))+
  xlab("Year")+
  ylab("Realised occupancy (psi.fs)")+
  annotate("text", x= 50, y= 0.9, parse=TRUE, label = paste("R[1]^{2}==" , round(cor(ret[[1]], simData$occSites), digits = 4)), size = 10)+
  annotate("text", x= 50, y= 0.8, parse=TRUE,label = paste("R[2]^{2}==" , round(cor(ret[[2]], simData$occSites[1:45]), digits = 4)), size = 10)+
  annotate("text", x= 50, y= 0.7, parse=TRUE,label = paste("R[3]^{2}==", round(cor(ret[[3]], simData$occSites), digits = 4)), size = 10)+
  annotate("text", x= 50, y= 0.6, parse=TRUE,label = paste("R[4]^{2}==" , round(cor(ret[[4]], simData$occSites), digits = 4)), size = 10)+
   theme(axis.title = element_text(size = 25),
         axis.text = element_text(size = 20),
         legend.title = element_text(size=30),
         legend.position = "bottom",
        legend.text = element_text(size=20),
        plot.title = element_text(size = 35))

ggsave(filename = "Figures/psifs.png",
       plot = occSites,
       width = 22,
       height = 10,
       dpi = 100)


## Extract time
ret <- lapply(allModels, function(x){
  x[[3]]
})

## Model Parameters and plot them
pars <- c("alphaP",           "alphaPSig"  ,      "alphaPhi" ,
 "alphaPsi"   ,      "alphaPsiSig" ,     "betaP" ,
 "betaPSig"  ,       "betaPhi"   ,       "betaPsi" ,
"betaPsiSig"    ,   "colonisationProb")

# Mean of pars
ret1 <- lapply(allModels, function(x){
  x[[2]]$all.chains[rownames(x[[2]]$all.chains) %in% pars, "Mean"]
})%>%
  do.call("cbind", .)%>%
  as.data.frame()%>%
  mutate(truth = c(4.99, 2,
                   2,
                   -0.33, 2,
                   0.227, 3,
                   2,
                   -1.64, 2,
                   0.05
  ),
  metric = "mean",
  parameters = pars)
colnames(ret1)[1:4] <- c('BMC',
                         'RMC', 'BUMC','AUMC')


# sd of pars
ret2 <- lapply(allModels, function(x){
  x[[2]]$all.chains[rownames(x[[2]]$all.chains) %in% pars, 3]
})%>%
  do.call("cbind", .)%>%
  as.data.frame()%>%
  mutate(truth = NA,
  metric = "sd",
  parameters = pars)
colnames(ret2)[1:4] <- c('BMC',
                         'RMC', 'BUMC','AUMC')

parEsts <- rbind(ret1,
                 ret2)

write.csv(parEsts, file = "Figures/occupancyParsEst.csv", row.names = FALSE)


## Effective sample size
ESSret <- lapply(allModels, function(x) {
  ret <- mcmcse::ess((as.matrix(x$samples))[,pars])%>%
    as.data.frame()
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(ESSret)[1:4] <- c('BMC',
                         'RMC', 'BUMC','AUMC')

ESSret1 <- ESSret%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  ggplot(., aes(x = parameters, y = value, col = variable))+
  geom_point(aes(shape = variable), size = 8,position = position_dodge(width = 0.9))+
  theme_classic()+
  scale_x_discrete(labels = c('alphaPSig' = expression(sigma[alpha^p]),
                              'betaPSig'= expression(sigma[beta^p]),
                              'alphaPsiSig'= expression(sigma[alpha^psi]),
                              'betaPsiSig'= expression(sigma[beta^psi]),
                              'alphaPhi'= expression(alpha^phi),
                              'betaPhi'= expression(beta^phi),
                              'alphaP'= expression(alpha^p),
                              'betaP'= expression(beta^p),
                              'alphaPsi'= expression(alpha^Psi),
                              'betaPsi'= expression(beta^Psi),
                              "colonisationProb"= expression(gamma)))+
  scale_color_manual(values = fill.colors)+
  xlab("")+
  ylab("ESS")+
  #labs(title = "A) M = 100")+
  labs(col = "Model", shape = "Model")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        plot.title = element_text(size = 35))


## Efficiency

EfficiencyRet <- lapply(allModels, function(x) {
  nDim <- length(x$timeRun)

  ess <- mcmcse::ess((as.matrix(x$samples))[,pars])%>%
    as.data.frame()
  timesRet <- if(nDim == 1){
    as.numeric(x$timeRun, units = "secs")
  }else{
    as.numeric(x$timeRun$all.chains, units = "secs")
  }

  ret <- ess/timesRet
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(EfficiencyRet)[1:4] <- c('BMC',
                                  'RMC', 'BUMC','AUMC')

EfficiencyRet1 <- EfficiencyRet%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  ggplot(., aes(x = parameters, y = value, col = variable))+
  geom_point(aes(shape = variable), size = 8,position = position_dodge(width = 0.9))+
  theme_classic()+
  scale_x_discrete(labels = c('alphaPSig' = expression(sigma[alpha^p]),
                              'betaPSig'= expression(sigma[beta^p]),
                              'alphaPsiSig'= expression(sigma[alpha^psi]),
                              'betaPsiSig'= expression(sigma[beta^psi]),
                              'alphaPhi'= expression(alpha^phi),
                              'betaPhi'= expression(beta^phi),
                              'alphaP'= expression(alpha^p),
                              'betaP'= expression(beta^p),
                              'alphaPsi'= expression(alpha^Psi),
                              'betaPsi'= expression(beta^Psi),
                              "colonisationProb"= expression(gamma)))+
  scale_color_manual(values = fill.colors)+
  xlab("")+
  ylab("Efficiency")+
  labs(col = "Model", shape = "Model")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        plot.title = element_text(size = 35))




## Put all plots together
fig <- ggpubr::ggarrange(ESSret1,
                         EfficiencyRet1,
                         ncol = 1, nrow = 2,
                         common.legend = TRUE,
                         legend = "top")

effESSPlot <- annotate_figure(fig,
                              bottom = text_grob("Model Parameters", size = 25))


ggsave(filename = "Figures/essEffPlotEx2.png",
       plot = effESSPlot,
       width = 18,
       height = 10,
       dpi = 100)

# Convergence of parameters
# Auxiliary PF
aPlots <- lapply(allModels, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% pars[1:4])%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots,
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D"))
aPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90),
                          bottom = text_grob("Iterations"))
ggsave(filename = "Figures/convergenceSecExA.png",
       plot = aPlots,
       width = 8,
       height = 8,
       units = "in")


aPlots <- lapply(allModels, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% pars[5:8])%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots,
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D"))
aPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90),
                          bottom = text_grob("Iterations"))
ggsave(filename = "Figures/convergenceSecExB.png",
       plot = aPlots,
       width = 8,
       height = 8,
       units = "in")

aPlots <- lapply(allModels, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% pars[9:11])%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots,
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D"))
aPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90),
                          bottom = text_grob("Iterations"))
ggsave(filename = "Figures/convergenceSecExC.png",
       plot = aPlots,
       width = 8,
       height = 8,
       units = "in")

timeRet <- lapply(allModels, function(x) {
  nDim <- length(x$timeRun)

  timesRet <- if(nDim == 1){
    as.numeric(x$timeRun, units = "mins")
  }else{
    as.numeric(x$timeRun$all.chains, units = "mins")
  }

  ret <- timesRet
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)
colnames(timeRet)[1:4] <- c('BMC',
                                  'RMC', 'BUMC','AUMC')

rm(list=ls())

###########################################
# Demographic SSM
##########################################

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)

# Fill colors for plotting
fill.colors <- c("BMC" = "#FF6600",
                 "ABSC" = "#FF3300",
                 "RMC" = "#00CCFF",
                 "BBSC" = "#E69F00",
                 "AUMC" = "#0033FF",
                 "BUMC" = "#3399FF",
                 "truth" = "black")

#bootstrap


load("demographicSSM/bootstrapPF/baselineModelMCMCResults.RData")
BBMC <- baselineModelMCMC
rm(baselineModelMCMC)

load("demographicSSM/bootstrapPF/reducedModelResults.RData")
BRMC <- example2ReducedModelTrue
rm(example2ReducedModelTrue)

load("demographicSSM/bootstrapPF/updatedModelResults.RData")
BUMC <- example2UpdatedModelTrue
rm(example2UpdatedModelTrue)

load("demographicSSM/auxilaryPF/updatedModelResults.RData")
AUMC <- example2UpdatedModelTrue
rm(example2UpdatedModelTrue)


# Put models together
allModels <- list(BBMC,
                  BRMC,
                  BUMC,
                  AUMC)

# Extract time
ret <- lapply(allModels, function(x){
  x[[3]]
})


# Extra growth rate and its standard deviation
ret <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("gammaX", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Median"]
})


retSD <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("gammaX", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, 3]
})


# Extract population index and its standard deviation
ret1 <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("popindex", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Median"]
})

ret1SD <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("popindex", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, 3]
})


# Plot population index
popnIndex <- data.frame(#BBSC = ret1[[1]],
  BMC = ret1[[1]],
  RMC = c(ret1[[2]], rep(NA, 5)),
  BUMC = c(ret1[[3]]),
  #ABSC = ret1[[5]],
  AUMC = c(ret1[[4]]),
  year = 1999:2016)%>%
  reshape2::melt(id.vars = c("year"),
                 value.name = "mean")

popnIndexSD <- data.frame(#BBSC = ret1[[1]],
  BMC = ret1SD[[1]],
  RMC = c(ret1SD[[2]], rep(NA, 5)),
  BUMC = c(ret1SD[[3]]),
  #ABSC = ret1[[5]],
  AUMC = c(ret1SD[[4]]),
  year = 1999:2016)%>%
  reshape2::melt(id.vars = c("year"))

popnIndex <- cbind(popnIndex, sd = popnIndexSD[,3])%>%
  ggplot(data = ., mapping = aes(x = year,
                                 y = mean, col = variable, group = variable))+
  geom_point(position = position_dodge(width = 0.7),size=4)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),linewidth = 1)+
  geom_line(position = position_dodge(width = 0.7),linewidth = 1)+
  annotate("text", x= 2005, y= 2500, parse=TRUE, label = paste("R[1]^{2}==" , round(cor(popnIndex[popnIndex$variable == "BUMC", 3][1:16],
                                                                                       popnIndex[popnIndex$variable == "BMC", 3][1:16]), digits = 4)), size = 10)+
  annotate("text", x= 2005, y= 2300, parse=TRUE,label = paste("R[2]^{2}==" , round(cor(popnIndex[popnIndex$variable == "AUMC", 3][1:16],
                                                                                      popnIndex[popnIndex$variable == "BMC", 3][1:16]), digits = 4)), size = 10)+
  #annotate("text", x= 2010, y= 0.4, parse=TRUE,label = paste("R[3]^{2}==", round(cor(allData[allData$Model == "RMC" & allData$variable =="psi.fs", 4][1:45],
   #                                                                                  allData[allData$Model == "BMC" & allData$variable =="psi.fs", 4][1:45]), digits = 4)), size = 10)+
  geom_vline(xintercept = 2011, linetype = "dashed", col = "red")+
  theme_classic()+
  xlab("Year")+
  ylab("Population Size")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))+
  scale_color_manual(values = fill.colors)

growthRate <- data.frame(#BBSC = ret[[1]],
  BMC = ret[[1]],
  RMC = c(ret[[2]], rep(NA, 5)),
  BUMC = c(ret[[3]]),
  #ABSC = ret[[5]],
  AUMC = c(ret[[4]]),
  year = 1999:2015)%>%
  reshape2::melt(id.vars = c("year"),
                 value.name = "mean")

growthRateSD <- data.frame(#BBSC = ret[[1]],
  BMC = retSD[[1]],
  RMC = c(retSD[[2]], rep(NA,5)),
  BUMC = c(retSD[[3]]),
  #ABSC = ret[[5]],
  AUMC = c(retSD[[4]]),
  year = 1999:2015)%>%
  reshape2::melt(id.vars = c("year"))

growthRate <- cbind(growthRate, sd = growthRateSD[,3])%>%
  ggplot(data = ., mapping = aes(x = year,
                                 y = mean, col = variable, group = variable))+
  geom_point(position = position_dodge(width = 0.7),size=4)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),linewidth = 1)+
  geom_line(position = position_dodge(width = 0.7),linewidth = 1)+
  theme_classic()+
  ylim(c(0.7, 1.6))+
  geom_hline(yintercept = 1)+
  annotate("text", x= 2005, y= 1.4, parse=TRUE, label = paste("R[1]^{2}==" , round(cor(growthRate[growthRate$variable == "BUMC", 3][1:15],
                                                                                       growthRate[growthRate$variable == "BMC", 3][1:15]), digits = 4)), size = 10)+
  annotate("text", x= 2005, y= 1.3, parse=TRUE,label = paste("R[2]^{2}==" , round(cor(growthRate[growthRate$variable == "AUMC", 3][1:15],
                                                                                      growthRate[growthRate$variable == "BMC", 3][1:15]), digits = 4)), size = 10)+
  scale_color_manual(values = fill.colors)+
  xlab("Year")+
  ylab("Growth rate")+
  geom_vline(xintercept = 2011, linetype = "dashed", col = "red")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

fig <- ggarrange(popnIndex,
                 growthRate,
                 ncol = 1,
                 nrow= 2,
                 common.legend = TRUE,
                 legend = "top",
                 labels = c("A)", "B)"),
                 font.label = list(size = 25))

ggsave(filename = "Figures/example2.png",
       plot = fig,
       width = 18,
       height = 10,
       dpi = 100)


pars <- c('mean.lambda',
          'beta.lam',
          'mean.gamma',
          'beta.gam',
          'mean.p',
          'beta.p',
          'logrho',
          'sd.rho',
          'sd.gam',
          'sd.lam'
)

pars <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[rownames(x[[2]]$all.chains) %in%pars]
  x[[2]]$all.chains[extractNames, 1]
})%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = c("BBMC",
                   "BRMC",
                   "BUMC",
                   "AUMC"))

write.csv(pars, file = "Figures/demograhicPars.csv", row.names = FALSE)

rm(list = ls())


#############
# Sparta
###############

# load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)

fill.colors <- c("BMC" = "#FF6600",
                 "ABSC" = "#FF3300",
                 "RMC" = "#00CCFF",
                 "BBSC" = "#E69F00",
                 "AUMC" = "#0033FF",
                 "BUMC" = "#3399FF",
                 "truth" = "black")

# load data
load("/Users/kwakupa/Dropbox/Data for PHD/particleFilters/spartaOccupancyModel/reducedModelResults.RData")
BRM <- example2ReducedModelTrue

load("/Users/kwakupa/Dropbox/Data for PHD/particleFilters/spartaOccupancyModel/updatedModelResultsBootstrap.RData")
BUMC <- example2UpdatedModelTrue

#load("/Users/kwakupa/Dropbox/Data for PHD/particleFilters/spartaOccupancyModel/updatedModelResultsAuxiliary.RData")
load("/Users/kwakupa/Dropbox/Data for PHD/particleFilters/spartaOccupancyModel/updatedModelResults.RData")
AUMC <- example2UpdatedModelTrue

rm(example2ReducedModelTrue, example2UpdatedModelTrue)

allModels <- list(BRM, BUMC, AUMC)


## Extract psi.fs
ret <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  rr <- x[[2]]$all.chains[extractNames, "Median"]
  if(length(rr) < 54){
    rr <- c(rr, rep(NA,5))
  }
  return(rr)
})%>%
  do.call("c", .)


retSD <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  rr <- x[[2]]$all.chains[extractNames, 3]
  if(length(rr) < 54){
    rr <- c(rr, rep(NA,5))
  }
  return(rr)
})%>%
  do.call("c", .)


load("spartaOccupancyModel/baselineModel.RData")

names <- rownames(baselineModel[[2]]$all.chains)[grepl("psi.fs", rownames(baselineModel[[2]]$all.chains))]
psi.fs <- baselineModel[[2]]$all.chains[names, "Mean"]

psi.fsSD <- baselineModel[[2]]$all.chains[names, 3]

# Extract a

retA <- lapply(seq_along(allModels), function(x){
  if(x == 1){
    vals <- c(allModels[[x]]$summary$all.chains[1:49, 1], rep(NA, 5))
  }else{
    vals <- allModels[[x]]$summary$all.chains[1:54, 1]
  }
})%>%
  do.call("c", .)

retAsd <- lapply(seq_along(allModels), function(x){
  if(x == 1){
    vals <- c(allModels[[x]]$summary$all.chains[1:49, 3], rep(NA,5))
  }else{
    vals <- allModels[[x]]$summary$all.chains[1:54, 3]
  }
})%>%
  do.call("c", .)

aOut <- baselineModel$summary$all.chains[1:54, 1]

aOutsd <- baselineModel$summary$all.chains[1:54, 1]
#aOut <- out$BUGSoutput$mean$a

year <- 1970:2023

extractedValuesA <- data.frame(a = c(retA, aOut),
                               psi.fs = c(ret, psi.fs),
                               #aSD = c(retAsd, aOutsd),
                               #psi.fsSD = c(retSD, psi.fsSD),
                               Model = rep(c("RMC","BUMC", "AUMC", "BMC"), each = 54),
                               year = rep(year, 4))%>%
  reshape2::melt(id.vars = c("Model", "year"),
                 value.name = "mean")


extractedValuesB <- data.frame(#a = c(retA, aOut),
  #psi.fs = c(ret, psi.fs),
  aSD = c(retAsd, aOutsd),
  psi.fsSD = c(retSD, psi.fsSD),
  Model = rep(c("RMC","BUMC", "AUMC", "BMC"), each = 54),
  year = rep(year, 4))%>%
  reshape2::melt(id.vars = c("Model", "year"),
                 value.name = "sd")

allData <- cbind(extractedValuesA, sd = extractedValuesB[,4 ])

#estimate the correlation
corA <- cor(allData[allData$Model == "BUMC" & allData$variable =="a", 4],
            allData[allData$Model == "BMC" & allData$variable =="a", 4])

corPsiFs <- cor(allData[allData$Model == "BUMC" & allData$variable =="psi.fs", 4],
            allData[allData$Model == "BMC" & allData$variable =="psi.fs", 4])

dat_text <- data.frame(
  label = c(label = paste("Cor =" , round(corA, digits = 2)), label = paste("Cor =" , round(corPsiFs, digits = 2))),
  variable   = c("a", "psi.fs"),
  Model = c("BUMC", "BMC"),
  x     = c(2018, 2018),
  y     = c(-8, 0.3)
)



extractedValuesA <- allData %>%
  filter(variable %in% c("a"))%>%
  filter(!Model %in% c("BUMC"))%>%
  ggplot(., mapping = aes(x = year, y = mean, col = Model))+
  geom_point(position = position_dodge(width = 0.7), size = 2)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),linewidth = 1)+
  geom_line(position = position_dodge(width = 0.7),linewidth = 1)+
  theme_classic()+
  ylab("Estimated value")+
  xlab("Year")+
  geom_vline(xintercept = 2018, linetype = "dashed")+
 # annotate("text", x= 2010, y= -8, parse=TRUE, label = paste("R[1]^{2}==" , round(cor(allData[allData$Model == "BUMC" & allData$variable =="a", 4],
  #                                                                                    allData[allData$Model == "BMC" & allData$variable =="a", 4]), digits = 4)), size = 10)+
  annotate("text", x= 2010, y= -7, parse=TRUE,label = paste("R[1]^{2}==" , round(cor(allData[allData$Model == "AUMC" & allData$variable =="a", 4],
                                                                                    allData[allData$Model == "BMC" & allData$variable =="a", 4]), digits = 4)), size = 10)+
  annotate("text", x= 2010, y= -10, parse=TRUE,label = paste("R[2]^{2}==", round(cor(allData[allData$Model == "RMC" & allData$variable =="a", 4][1:45],
                                                                                   allData[allData$Model == "BMC" & allData$variable =="a", 4][1:45]), digits = 4)), size = 10)+
  #annotate("text", x= 50, y= 0.6, parse=TRUE,label = paste("R[4]^{2}==" , round(cor(ret[[4]], simData$occSites), digits = 4)), size = 10)+
  #annotate("text", x= 2018, y= -8, label = paste("Cor =" , round(corA, digits = 2)))+
  #annotate("text", x= 2018, y= 0.3, label = paste("Cor =" , round(corPsiFs, digits = 2)))+
  facet_wrap( ~variable, ncol = 1, nrow = 2, scales = "free_y")+
  theme(legend.position = "top")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))+
  scale_color_manual(values = fill.colors)

extractedValuesPsiFs <- allData %>%
  filter(variable %in% c("psi.fs"))%>%
  filter(!Model %in% c("BUMC"))%>%
  ggplot(., mapping = aes(x = year, y = mean, col = Model))+
  geom_point(position = position_dodge(width = 0.7), size = 2)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),linewidth = 1)+
  geom_line(position = position_dodge(width = 0.7),linewidth = 1)+
  theme_classic()+
  ylab("Estimated value")+
  xlab("Year")+
  geom_vline(xintercept = 2018, linetype = "dashed")+
  annotate("text", x= 2010, y= 0.2, parse=TRUE,label = paste("R[1]^{2}==" , round(cor(allData[allData$Model == "AUMC" & allData$variable =="psi.fs", 4],
                                                                                     allData[allData$Model == "BMC" & allData$variable =="psi.fs", 4]), digits = 4)), size = 10)+
  annotate("text", x= 2010, y= 0.35, parse=TRUE,label = paste("R[2]^{2}==", round(cor(allData[allData$Model == "RMC" & allData$variable =="psi.fs", 4][1:45],
                                                                                     allData[allData$Model == "BMC" & allData$variable =="psi.fs", 4][1:45]), digits = 4)), size = 10)+

  facet_wrap( ~variable, ncol = 1, nrow = 2, scales = "free_y")+
  theme(legend.position = "top")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))+
  scale_color_manual(values = fill.colors)

extractedValues <- ggarrange(extractedValuesA,
                             extractedValuesPsiFs,
                             ncol=1)

ggsave(filename = "Figures/spartaPsifs.png",
       plot = extractedValues,
       width = 18,
       height = 10,
       dpi = 100)
