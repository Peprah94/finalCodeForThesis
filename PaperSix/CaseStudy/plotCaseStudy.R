source("idmFunctions/functionForEstimation.R")
library(readr)
library(dplyr)
load("CaseStudy/data_idm.RData")

#Extracting values from formatted data for Table 1 in main paper
lapply(simulations_all, function(x){
  sumIsNA <- function(z){
    sum(!is.na(c(z)))/nrow(z)
  }
  reportingRate <- mean(apply(x$mat.species, c(3), function(x) mean(x> 0, na.rm = TRUE)), na.rm = TRUE)
  reportingSD <- sd(apply(x$mat.species, c(3), function(x) mean(x> 0, na.rm = TRUE)), na.rm = TRUE)
  NaiveOccupancy <-  mean(apply(x$mat.species, c(1,2), function(x) mean(x> 0, na.rm = TRUE)), na.rm = TRUE)
  NaiveOccupancySD <-  sd(apply(x$mat.species, c(1,2), function(x) mean(x> 0, na.rm = TRUE)), na.rm = TRUE)
   #sum(x$mat.species > 0 , na.rm = TRUE)/length(!is.na(x$mat.species))
  missSitesSCounts <- sum(x$mat.genus, na.rm = TRUE)/length(!is.na(x$mat.genus))
  missSitesSpeciesSD <- sd(x$mat.species > 0, na.rm = TRUE)
  missSitesCountsSD <- sd(x$mat.genus, na.rm = TRUE)
  ret <- c(reportingRate, NaiveOccupancy, missSitesSCounts,
           reportingSD, NaiveOccupancySD,missSitesCountsSD)
  return(ret)
})


# Define parameters
models <- c("IDMSH" , "SOM", "IDMCO", "GCMSH", "GCMCO")
insectGrps <- c("Bumblebees", "Hoverflies", "Solitary bees")
nGroups = length(insectGrps)
nModels = length(models)

#########################
# Plot estimates from MCMC
##########################

#load data from Figures folder
parsEstimates <- read_csv("Figures/sigmaEstimates.csv")
shannonEstimates <- read_csv("Figures/shannonEstimates.csv")
rHatEstimates <- read_csv("Figures/rHatEstimates.csv")

# plot convergence plots
rHatSigma <- rHatEstimates%>%
 # filter(!grepl("Cov|shan|eta.lam",Parameter))%>%
  filter(grepl("sigma|rho|mu",Parameter))%>%
  mutate(method = rep(rep(models, each = 150/(nGroups*nModels)), nGroups),
         insects = rep(insectGrps,
                       each= 150/nGroups))%>%
  ggplot(aes(x = Rhat, y = Parameter, col = method))+
  geom_point(position=position_dodge(0.5))+
  geom_vline(xintercept = 1.1, linetype = "dashed")+
  facet_grid( ~ insects)+
  theme_classic()+
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.placement = 'outside',
        strip.background = element_rect(colour = "black", fill = "white")#,
       # text = element_text(size = 20)
        )

ggsave(filename = "Figures/rHat.png",
       rHatSigma,
       width = 10,
       height = 7,
       dpi = 700)


# Plot mu estimates
muEstimates <- parsEstimates%>%
  #as.data.frame()%>%
  filter(!grepl("rho|sigma",variable))%>%
  dplyr::mutate(value = if_else(model %in% c("IGCO") & variable %in% c("muSpecies"), NA, value))%>%
  mutate(value = if_else(model %in% c("IGSH", "IDMSH", "SOM") & variable %in% c("mulambda"), NA, value)
         )%>%
  ggplot( aes(x=insects, y=value, col=as.factor(model),  shape = as.factor(model))) +
  geom_point(position=position_dodge(0.5))+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0.5,
                position=position_dodge(0.5))+
  facet_grid( ~ variable)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Insect group")+
  ylab(expression(paste("Posterior mean" %+-% "SE")))+
  theme_classic()+
  coord_trans()+
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.placement = 'outside',
        strip.background = element_rect(colour = "black", fill = "white"))+
  scale_color_discrete(labels=c('IDMCO', 'IDMSH', 'GCMCO', 'GCMSH', 'SOM'))+
  scale_shape_discrete(labels=c('IDMCO', 'IDMSH', 'GCMCO', 'GCMSH','SOM'))
ggsave(filename = "Figures/muEstimates.png",
       muEstimates,
       width = 10,
       height = 7,
       dpi = 700)
  #sigma estimates
  #remove "sigmaSpecies" after new results
 sigmaEstimates <-  parsEstimates%>%
    filter(!grepl("rho|mu",variable))%>%
   dplyr::mutate(value = if_else(model %in% c("IGCO") & variable %in% c("muSpecies"), NA, value))%>%
   mutate(value = if_else(model %in% c("IGSH", "IDMSH", "SOM") & variable %in% c("mulambda"), NA, value)
   )%>%
    #mutate(value = ifelse((variable %in% c("sigmaAlphaSites", "sigmaAlphaSpecies", "sigmaSpecies") & method == "Group only"), NA, value))%>%
    ggplot( aes(x=insects, y=value, col=as.factor(model), shape=as.factor(model))) +
    geom_point(position=position_dodge(0.5))+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.5))+
    facet_grid( ~ variable, scales = "free_y")+
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Insect group")+
    ylab(expression(paste("Posterior mean" %+-% "SE")))+
  theme_classic()+
    coord_trans()+
    theme(legend.position = "bottom", legend.title = element_blank(),
          strip.placement = 'outside',
          strip.background = element_rect(colour = "black", fill = "white"))#+
    #scale_color_discrete(labels=c('IDMCO', 'IDMSH', 'SOM'))+
    #scale_shape_discrete(labels=c('IDMCO', 'IDMSH', 'SOM'))
 ggsave(filename = "Figures/sigmaEstimates.png",
        sigmaEstimates,
        width = 10,
        height = 7,
        dpi = 700)


# rho estimates
 rhoEstimates <-  parsEstimates%>%
    filter(!grepl("sigma|mu",variable))%>%
    filter(model != "SOM")%>%
    #mutate(value = ifelse(method == "Species Only", NA, value))%>%
    ggplot( aes(x=insects, y=value, col=as.factor(model), shape = as.factor(model))) +
    geom_point(position=position_dodge(0.5))+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.5))+
    #facet_grid(. ~ model)+
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Insect group")+
    ylab(expression(paste("Posterior mean" %+-% "SE")))+
    theme_classic()+
    coord_trans()+
    theme(legend.position = "bottom", legend.title = element_blank(),
          strip.placement = 'outside',
          strip.background = element_rect(colour = "black", fill = "white"))#+
   # scale_color_discrete(labels=c('GCM', 'IDM', 'SOM'))+
    #scale_shape_discrete(labels=c('GCM', 'IDM', 'SOM'))

 ggsave(filename = "Figures/rhoEstimates.png",
        rhoEstimates,
        width = 10,
        height = 7,
        dpi = 700)

#shannon Index

  averageShannon <- shannonEstimates%>%
    group_by(insects, model) %>%
    summarise(average = mean(mean),
              sd = sd(mean))%>%
    ungroup()%>%
    ggplot( aes(x=insects, y=average, col=as.factor(model))) +
    geom_point(position=position_dodge(0.5))+
    geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                  position=position_dodge(0.5))+
    #facet_grid(~model)+
    scale_x_discrete(guide = guide_axis(angle = 0)) +
    xlab("Insect group")+
    ylab(expression(paste("Posterior mean" %+-% "SE")))+
    theme_article()+
    coord_trans()+
    scale_color_discrete(labels=c('IDMCO', 'IDMSH', 'GCMCO', 'GCMSH','SOM'), name = "Models")+
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          strip.placement = 'outside',
          strip.background = element_rect(colour = "black", fill = "white"))



  #save plot
  ggsave("Figures/averageShannonPlot.svg",
         plot = averageShannon,
         dpi = 320,
         height = 8,
         width = 10,
         units = "cm")

  # shanEstimates0 <- read_csv("simEstimates/shanSimulationsEstimates0.csv")
  # shanEstimates1 <- read_csv("simEstimates/shanSimulationsEstimates1.csv")
  # shanEstimates2 <- read_csv("simEstimates/shanSimulationsEstimates2.csv")
  #
  # shanEstimates <- rbind(shanEstimates0,
  #                        shanEstimates1,
  #                        shanEstimates2)
  #
  # shanEst <- shanEstimates%>%
  #   ggplot()+
  #   geom_boxplot(aes(x = model, y = V1, fill = truth))+
  #   ylim(c(0, 1.5))+
  #   # facet_grid(~truth)+
  #   ylab("MSE(H')")+
  #   theme_classic()+
  #   theme(legend.title = element_blank(),
  #         legend.position = "bottom",
  #         legend.text = element_text(size = 8))
  #
  #
  # gg1 <- ggpubr::ggarrange(averageShannon, shanEst,
  #                          labels = c("A)", "B)"),
  #                          nrow = 2, ncol = 1,legend = "top")
  #
  # ggsave("Figures/simHills.svg",
  #        plot = gg1,
  #        dpi = 720,
  #        height = 14,
  #        width = 10,
  #        units = "cm")

  # Plot the Shannon Indices on the map
  # map of sites (represented at 10Km because too small to see at 1Km)
  library(BRCmap)
  coords <- gr_let2num(shannonEstimates$coords)
  shannonEstimates <- cbind(coords, shannonEstimates)
  UK_countries_lowres@data$ID = rownames(UK_countries_lowres@data)
  UK_countries_lowres_df <- ggplot2::fortify(UK_countries_lowres)

  #plot Shannon map for each group
  shannonMap <- lapply(as.list(c("Bumblebees", "Hoverflies" , "Solitary bees")), function(x){
    shannonEstimates%>%
      filter(insects == x)%>%
       mutate(model = ifelse(model == "IGCO", "GCMCO",
                             ifelse(model == "IGSH", "GCMSH", model)))%>%
   ggplot() +
      geom_polygon(data = UK_countries_lowres_df,
                   aes(x = long, y = lat, group = group),
                   fill = "gray93", colour = "black") +
      geom_point(aes(x = EASTING, y = NORTHING,
                     fill = mean), shape = 22, size = 1,
                 color = "gray") +
      facet_grid(~model) +
      cowplot::theme_map(font_size = 10) +
      theme(legend.position = "right",
            #legend.title = element_blank(),
            strip.placement = 'outside',
            strip.background = element_rect(colour = "black", fill = "white"),
            text = element_text(size = 10))+
      #coord_equal() +
      scale_fill_viridis_c(name = "H'")
  })


  shannonMapRet <-  ggarrange(plots = shannonMap,
            nrow = 3,
            legend = "right",
           # labels = c("a) Bumblebees", "b) Hoverflies", "c) Solitary bees"),
           labels = c("A)", "B)", "C)"),
            label.args = list(gp = grid::gpar(font = 2, cex =
                                                1)))

  ggsave("Figures/shannonMapRet.svg",
         plot = shannonMapRet,
         dpi = 720,
         height = 18,
         width = 16,
         units = "cm")


  shannonMapSD <- lapply(as.list(c("Bumblebees", "Hoverflies" , "Solitary bees")), function(x){
    shannonEstimates%>%
      filter(insects == x)%>%
      # mutate(method = ifelse(method == "Group only", "GCM",
      #                        ifelse(method == "Species Only", "SOM", method)))%>%
      ggplot() +
      geom_polygon(data = UK_countries_lowres_df,
                   aes(x = long, y = lat, group = group),
                   fill = "gray93", colour = "black") +
      geom_point(aes(x = EASTING, y = NORTHING,
                     fill = sd), shape = 22, size = 2,
                 color = "gray") +
      facet_grid(~ model) +
      cowplot::theme_map(font_size = 10) +
      theme(legend.position = "right",
            #legend.title = element_tex(),
            strip.placement = 'outside',
            strip.background = element_rect(colour = "black", fill = "white"),
            text = element_text(size = 10))+
      coord_equal() +
      scale_fill_viridis_c(name = "SD(H')")
  })

  shannonMapSDRet <-  ggarrange(plots = shannonMapSD,
                              nrow = 3,
                              legend = "right",
                              # labels = c("a) Bumblebees", "b) Hoverflies", "c) Solitary bees"),
                              labels = c("A)", "B)", "C)"),
                              label.args = list(gp = grid::gpar(font = 2, cex =
                                                                  1)))

  ggsave("Figures/shannonMapSDRet.png",
         plot = shannonMapSDRet,
         dpi = 320,
         height = 12,
         width = 10,
         units = "cm")





  ###########################
  # Plot posterior predictive checks
  ###########################

  #load data
  load("CaseStudy/bumblebees/posteriorPredCheckSmall.RData")
  bb <- ret
  load("CaseStudy/hoverflies/posteriorPredCheckSmall.RData")
  hv <- ret
  load("CaseStudy/solitarybees/posteriorPredCheckSmall.RData")
  sb <- ret

  allResults <- c(bb,
                  hv,
                  sb)

  # return estimated p-values
 estPvalues <-  lapply(allResults, function(x){
    colMeans(x[,1:4])
  })%>%
    do.call("rbind",.)%>%
    as.data.frame("")%>%
    mutate(model = rep(c("IDMSH", "SOM", "IDMCO", "GCMSH", "GCMCO"), 3),
           insectGrps = rep(insectGrps, each = 5))
colnames(estPvalues)[1:4] <- c("PanTrapDeviance", "FITcntDev", "PanTrapMean", "FITcntMean")


## Histogram of values

# histPlots <- lapply(allResults, function(x){
#   #convert data to data frame
#   x <- x%>%
#     as.data.frame()
#
#   #pantrap data
#   ## Extract x-intercept, p-value and points to annotate figure
#   x_intercept <-x$V4[1]
#   x1 <- max(x$V3) - 0.005
#   y1 <- 5
#   bpv <- mean(x$V1)
#   if(is.nan(bpv)) bpv <- "NA"
#   pvalText <- paste("Bpv = ", bpv)
#
#   ##plot data
#   pantrap <- x %>%
#     ggplot()+
#     geom_histogram(mapping = aes(x = V3))+
#     ylab("Counts")+
#     xlab("Average occupied sites")+
#     theme_classic()+
#     geom_vline(xintercept = x_intercept, col = "red")+
#     annotate(geom = 'text', x = x1, y = y1,
#              label = pvalText)
#
#   #FIT count data
#   ## Extract x-intercept, p-value and points to annotate figure
#   x_intercept <-x$V8[1]
#   x1 <- max(x$V7) - 0.005
#   y1 <- 5
#   bpv <- mean(x$V9)
#   pvalText <- paste("Bpv = ", bpv)
#
#   ##plot data
#   fitCnts <- x %>%
#     ggplot()+
#     geom_histogram(mapping = aes(x = V7))+
#     ylab("Counts")+
#     xlab("Average FIT counts")+
#     theme_classic()+
#     geom_vline(xintercept = x_intercept, col = "red")+
#     annotate(geom = 'text', x = x1, y = y1,
#              label = pvalText)
#
#   ggpubr::ggarrange(pantrap, fitCnts, nrow = 1)
# })
#
# ggpubr::ggarrange( plotlist = histPlots, ncol = 4, nrow = 5,
#                    hjust = -2,
#                    labels = c("BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH",
#                               "BB - IDMSH"),
#                    label.args = list(gp = grid::gpar(font = 2, cex =
#                                                        1)))
#

# Check the species for each Bumblebees species
meanDevbb <- lapply(bb[1:3], function(x){
  dat <- colMeans(x[, 13:29])
  names <- as.factor(1:17)
  allData <- data.frame(dat = dat,
                        names = names)%>%
    mutate(condition = ifelse(dat > 0.4 & dat < 0.6, 1,0 ))
} )%>%
  do.call("rbind", .)%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO"), each = 17),
         statistic = "mean",
         insect = "bumblebees")



devianceDevbb <- lapply(bb[1:3], function(x){
  dat <- colMeans(x[, 30:46])
  names <- as.factor(1:17)
  allData <- data.frame(dat = dat,
                        names = names)%>%
    mutate(condition = ifelse(dat > 0.4 & dat < 0.6, 1,0 ))
} )%>%
  do.call("rbind", .)%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO"), each = 17),
         statistic = "deviance",
         insect = "bumblebees")

# devianceDevsb%>%
#   group_by(model) %>%
#   summarise(average = mean(condition),
#             sd = sd(condition))%>%
#   ungroup()
# allData <- rbind(meanDev,
#                  devianceDev)%>%
#   ggplot(., aes(x = names, y = dat, col = as.factor(condition)))+
#     geom_point()+
#     geom_hline(yintercept = 0.5, col = "red")+
#     geom_hline(yintercept = 0.6, linetype="dashed", col = "red")+
#     geom_hline(yintercept = 0.4, linetype="dashed", col = "red")+
#   theme_article()+
#     scale_color_manual(values = c("#E69F00", "#999999"), name = "Fit", labels = c("Bad", "Good"))+
#   facet_grid(statistic~model)+
#   xlab("species")+
#   ylab("Posterior predictive value")+
#   #scale_color_discrete(name = "Fit", labels = c("Bad", "Good"))+
#     theme(legend.position = "bottom")#


# Check the species for each Bumblebees species
meanDevhv <- lapply(hv[1:3], function(x){
  dat <- colMeans(x[, 13:91])
  names <- as.factor(1:79)
  allData <- data.frame(dat = dat,
                        names = names)%>%
    mutate(condition = ifelse(dat > 0.4 & dat < 0.6, 1,0 ))
} )%>%
  do.call("rbind", .)%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO"), each = 79),
         statistic = "mean",
         insect = "hoverflies")

devianceDevhv <- lapply(hv[1:3], function(x){
  dat <- colMeans(x[, 92:170])
  names <- as.factor(1:79)
  allData <- data.frame(dat = dat,
                        names = names)%>%
    mutate(condition = ifelse(dat > 0.4 & dat < 0.6, 1,0 ))
} )%>%
  do.call("rbind", .)%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO"), each = 79),
         statistic = "deviance",
         insect = "hoverflies")

# allData <- rbind(meanDev,
#                  devianceDev)%>%
#   ggplot(., aes(x = names, y = dat, col = as.factor(condition)))+
#   geom_point()+
#   geom_hline(yintercept = 0.5, col = "red")+
#   geom_hline(yintercept = 0.6, linetype="dashed", col = "red")+
#   geom_hline(yintercept = 0.4, linetype="dashed", col = "red")+
#   theme_article()+
#   scale_color_manual(values = c("#E69F00", "#999999"), name = "Fit", labels = c("Bad", "Good"))+
#   facet_grid(statistic~model)+
#   xlab("species")+
#   ylab("Posterior predictive value")+
#   #scale_color_discrete(name = "Fit", labels = c("Bad", "Good"))+
#   theme(legend.position = "bottom")#


# Check the species for each solitarybees species
meanDevsb <- lapply(sb[1:3], function(x){
  dat <- colMeans(x[, 13:82])
  names <- as.factor(1:70)
  allData <- data.frame(dat = dat,
                        names = names)%>%
    mutate(condition = ifelse(dat > 0.4 & dat < 0.6, 1,0 ))
} )%>%
  do.call("rbind", .)%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO"), each = 70),
         statistic = "mean",
         insect = "solitarybees")

devianceDevsb <- lapply(sb[1:3], function(x){
  dat <- colMeans(x[, 83:152])
  names <- as.factor(1:70)
  allData <- data.frame(dat = dat,
                        names = names)%>%
    mutate(condition = ifelse(dat > 0.4 & dat < 0.6, 1,0 ))
} )%>%
  do.call("rbind", .)%>%
  mutate(model = rep(c("IDMSH" , "SOM", "IDMCO"), each = 70),
         statistic = "deviance",
         insect = "solitarybees")

devianceDevsb%>%
  group_by(model) %>%
  summarise(average = mean(condition),
            sd = sd(condition))%>%
  ungroup()

allData <- rbind(devianceDevbb,
                 devianceDevhv,
                 devianceDevsb)%>%
  ggplot(., aes(x = names, y = dat, col = as.factor(condition)))+
  geom_point()+
  geom_hline(yintercept = 0.5, col = "red")+
  geom_hline(yintercept = 0.6, linetype="dashed", col = "red")+
  geom_hline(yintercept = 0.4, linetype="dashed", col = "red")+
  xlab("")+
  theme_article()+
  scale_color_manual(values = c("#E69F00", "#999999"), name = "Fit", labels = c("Bad", "Good"))+
  facet_grid(model~insect, scales = "free")+
  xlab("species")+
  ylab("Posterior predictive value")+
  #scale_color_discrete(name = "Fit", labels = c("Bad", "Good"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, size = 4))#

ggsave("Figures/posteriorPredCheck.png",
       plot = allData,
       dpi = 320,
       height = 10,
       width = 16,
       units = "cm")
