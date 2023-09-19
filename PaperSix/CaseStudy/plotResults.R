setwd("/Volumes/kwakupa-1/IDM_new/Data")
load("/Volumes/kwakupa-1/IDM_new/Data/data_idm.RData")
source("functionForEstimation.R")

#set theme for plots
theme_set(theme_classic())


# Function that takes the name of the folder and 
#generates all the figures needed for the paper

DataFormatting <- function(name){
  
#loading data from folder
load(paste(name,"estimate_data_inter_bum1.RData", sep="/"))
bb_inter <- estimates
load(paste(name,"estimate_data_species_bum1.RData", sep="/"))
bb_spe <- estimates
load(paste(name,"estimate_data_genus_bum1.RData", sep="/"))
bb_gen <- estimates
load(paste(name,"estimate_data_inter_hov1.RData", sep="/"))
hov_inter <- estimates
load(paste(name,"estimate_data_species_hov1.RData", sep="/"))
hov_spe <- estimates
load(paste(name,"estimate_data_genus_hov1.RData", sep="/"))
hov_gen <- estimates
load(paste(name,"estimate_data_inter_sol1.RData", sep="/"))
sol_inter <- estimates
load(paste(name,"estimate_data_species_sol1.RData", sep="/"))
sol_spe <- estimates
load(paste(name,"estimate_data_genus_sol1.RData", sep="/"))
#load("shared/estimate_data_genus_sol1.RData")
sol_gen <- estimates

#Put all data together
all_data <- list(bb_inter , bb_spe , bb_gen ,
              hov_inter , hov_spe , hov_gen,
              sol_inter ,  sol_spe , sol_gen )

# Retrieving names of species
speciesNamesList <- c(rep(list(simulations_all[[1]]$species_names),3),
                rep(list(simulations_all[[2]]$species_names),3),
                rep(list(simulations_all[[3]]$species_names),3)
)


#retrieving constants and adding it to the MCMCoutput data
simData = rep(simulations_all, each = 3) #Data used for MCMC

allData <- lapply(as.list(seq(1,9)), function(x){
  
  # Extracting dimensions of data parameters
  data <- simData[[x]]
  dim_data <- dim(data[[1]]) 
  data_dim <- dim_data[2] 
  
  # Constants used for the MCMC
  const <- list(n.sites = dim_data[1], 
                n.species= (dim_data[2]), 
                n.years = dim_data[4],
                n.replicates = 5,
                n.visit=dim_data[3],
                mu.eta = rep(0, dim_data[2]),
                nu = 3,
                NLatent=ceiling(data_dim/10),
                a1 = 50,
                a2= 50,
                a.sigma = 2,
                b.sigma = 1
  )

  retData <- c(const, all_data[[x]])
  return(retData)
  
})




 #Abbreviate the names of the species
abbreSpeciesNamesList <- lapply(speciesNamesList, function(x){
  ret <- abbreviate(x, minlength = 8)
  sub("(\\w+\\s+\\w+).*", "\\1", ret)
})
 

##################
#   Species Interraction plot
##############
 correlationMatrices <- pmap(list(allData, 
                                  abbreSpeciesNamesList), 
                             extractCovariance)

title_text <- c("a) bumblebees \n IDM", 
                "b) bumblebees \n Speices only", 
                "c) bumblebees \n Group only",
                   "d) hoverflies \n IDM",
                "e) hoverflies \n Speices only", 
                "f) hoverflies \n Group only",
                   "g) solitarybees \n IDM", 
                "h) solitarybees \n Speices only", 
                "i) solitarybees \n Group only")
 
#Plot the correlation matrices
  ggCorrplots <-  lapply(seq_along(correlationMatrices), function(x){ 
     ggcorrplot(correlationMatrices[[x]],
              hc.order = TRUE,
              type = "upper",
              tl.cex = 2,
              tl.srt = 90,
              lab_size = 2,
              title = title_text[[x]])+
      theme(plot.title = element_text(size = 9) )
              }) 

  # Arrange the plots
   gg2 <- ggpubr::ggarrange(plotlist = ggCorrplots,
                            ncol = 3, nrow = 3, common.legend = TRUE, 
                            legend = "right",
                            heights = c(1.5,2,2)) 
  
   figure <- ggpubr::annotate_figure(gg2, left = "Species list", bottom = "Species list")
   
   #save plot
   ggsave("Figures/species_association.png", plot = figure, dpi = 320, 
          height = 18, width = 18, units = "cm")

      
 ##########################################
 # extracting the sd of the site effects
 ############################################
 sigmaEstimates <- lapply(allData, extractSigma)%>%
   do.call("rbind", .)%>%
   as.data.frame()%>%
   mutate(method = rep(c("IDM" , "Species only", "Group Only"), 3),
          insects = rep(c("Bumblebees", "Solitary bees", "Hoverflies"),
                        each= 3))

colnames(sigmaEstimates) <- c("Site detection", "Site Observed", "Methods", "Insect Group")
write.csv(sigmaEstimates,"Figures/sigmaEstimates.csv", row.names = FALSE)

# Plot the standard deviations
plotSigmas <- melt(sigmaEstimates, id.vars = c("Methods", "Insect Group"))%>%
  ggplot(aes(x= Methods, y = value, group = `Insect Group`, col= `Insect Group`))+
  geom_point()+
  geom_line(linetype="dashed")+
  theme_bw()+
  ylab("Variance")+
  facet_wrap(~variable, nrow=2)

#Save plot
ggsave(filename = "Figures/sdEstimates.png",plotSigmas)


#########################
# Beta0 estimates
########################
#Extreact bete0 values
beta0Matrices <- pmap(list(allData, abbreSpeciesNamesList), extractBeta0)
beta0SDMatrices <- pmap(list(allData, abbreSpeciesNamesList), extractBeta0SD)
#Title for their plot
titleTextBeta0 <- c("a) bumblebees \n IDM", "b) bumblebees \n Speices only", "c) bumblebees \n Group only",
                "d) hoverflies \n IDM", "e) hoverflies \n Speices only", "f) hoverflies \n Group only",
                "g) solitarybees \n IDM", "h) solitarybees \n Speices only", "i) solitarybees \n Group only")

#Plot beta0
ggBeta0 <-  lapply(seq_along(beta0Matrices), function(x){ 
beta0Matrices[[x]]%>%
    arrange(mean)%>%
    mutate(id = row_number(),
           species = factor(species, levels = unique(species)))%>%
  ggplot()+
    geom_point(aes(y = species, x = mean))+
    #scale_y_discrete(breaks = x$id, labels = x$species)+
    geom_segment(aes(y = species,yend = species,
                     x = low, xend = upper))+
    labs(title = titleTextBeta0[[x]])+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(plot.title = element_text(size = 7),
          axis.text.x = element_text(size = 8), 
          axis.text.y =element_text(size = 3) )
}) 

#Plot for SD of beta0
ggbeta0SD <- lapply(seq_along(beta0SDMatrices), function(x){ 
  beta0SDMatrices[[x]]%>%
    arrange(precision)%>%
    mutate(id = row_number(),
           precision = 1/precision,
           species = factor(species, levels = unique(species)))%>%
    #arrange(precision)%>%
    ggplot()+
    geom_point(aes(y = species, x = precision))+
    labs(title = titleTextBeta0[[x]])+
    xlab("")+
    ylab("")+
    theme_bw()+
    #scale_x_continuous(trans='log')+
    theme(plot.title = element_text(size = 7),
          axis.text.x = element_text(size = 8), 
          axis.text.y =element_text(size = 3) )
}) 

#Put all the plots together
gg1 <- ggpubr::ggarrange(plotlist = ggBeta0,
                         ncol = 3, nrow = 3, 
                common.legend = TRUE, 
                legend = "right",
                heights = c(2,3, 3)) 

gg1SD <- ggpubr::ggarrange(plotlist =ggbeta0SD,
                           ncol = 3, nrow = 3, 
                           common.legend = TRUE, 
                           legend = "right",
                           heights = c(2,3, 3)) 

figureBeta0 <- ggpubr::annotate_figure(gg1, 
                                       left = "Species", 
                                       bottom = "Beta0 estimate")
figureBeta0SD <- ggpubr::annotate_figure(gg1SD, 
                                       left = "Species", 
                                       bottom = "Beta0 estimate")
#save plot
ggsave("Figures/beta0.png", 
       plot = figureBeta0, 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")

ggsave("Figures/beta0SD.png", 
       plot = figureBeta0SD, 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")

#########################
# alpha0
###################

#Estreact alpha0 estimates
alpha0Matrices <- pmap(list(allData, abbreSpeciesNamesList), 
                       extractAlpha0)

#plot data
ggAlpha0 <-  lapply(seq_along(alpha0Matrices), function(x){ 
  alpha0Matrices[[x]]%>%
    arrange(mean)%>%
    mutate(id = row_number(),
           species = factor(species, levels = unique(species)))%>%
    ggplot()+
    geom_point(aes(y = species, x = mean))+
    #scale_y_discrete(breaks = x$id, labels = x$species)+
    geom_segment(aes(y = species,yend = species,
                     x = low, xend = upper))+
    labs(title = titleTextBeta0[[x]])+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(plot.title = element_text(size = 7),
          axis.text.x = element_text(size = 8), 
          axis.text.y =element_text(size = 3) )
}) 

#put plots together
gg3 <- ggpubr::ggarrange(plotlist = ggAlpha0,
                         ncol = 3, nrow = 3, 
                         common.legend = TRUE, 
                         legend = "right",
                         heights = c(2,2, 3)) 

figureAlpha0 <- ggpubr::annotate_figure(gg3, 
                                        left = "Species", 
                                        bottom = "Intercept estimate")

ggsave("Figures/Alpha0.png", 
       plot = figureAlpha0, 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")


#############
#   Intercept Lambda
#############3

interceptLambda <- lapply(allData, function(x){
                            extractInterceptLambda(x)})%>%
                            do.call("rbind", .)

interceptLambdaSD <- lapply(allData, function(x){
  extractInterceptLambdaSD(x)})%>%
  do.call("rbind", .)

#take out data that belongs to the species 
idx <- c(2,5,8)

interceptLambda <- interceptLambda[-idx,]%>%
                    data.frame()%>%
                    mutate(method = rep(c("IDM", "Group only"), 3),
                      group = rep(c("bumblebees", "hoverflies", "solitarybees"),
                                  each = 2))

interceptLambdaSD <- interceptLambdaSD[-idx,]%>%
  data.frame()%>%
  mutate(method = rep(c("IDM", "Group only"), 3),
         group = rep(c("bumblebees", "hoverflies", "solitarybees"),
                     each = 2))

colnames(interceptLambda) <- c("mean", "low", "high", "method", "group")
colnames(interceptLambdaSD) <- c("precision", "method", "group")

#plot data
ggInterceptLambda <-  interceptLambda%>%
  mutate(group_method = c("bumblebees-IDM","bumblebees-Group Only" ,
                          "hoverflies-IDM","hoverflies-Group Only",
                          "solitarybees-IDM","solitarybees-Group Only"))%>%
  arrange(mean)%>%
  mutate(id = row_number(),
         species = factor(group_method, levels = unique(group_method)))%>%
  ggplot()+
  geom_point(aes(y = species, x = mean))+
  geom_segment(aes(y = species,yend = species,
                   x = low, xend = high))+
  xlab("gamma estimate")+
  ylab("Method & Insect Group")+
  theme_bw()+
  theme(plot.title = element_text(size = 7),
        axis.text.x = element_text(size = 8), 
        axis.text.y =element_text(size = 5) )

#Plot SD of intercept lambda
ggInterceptLambdaSD <-  interceptLambdaSD%>%
  mutate(group_method = c("bumblebees-IDM","bumblebees-Group Only" ,
                          "hoverflies-IDM","hoverflies-Group Only",
                          "solitarybees-IDM","solitarybees-Group Only"))%>%
  arrange(precision)%>%
  mutate(id = row_number(),
         species = factor(group_method, levels = unique(group_method)))%>%
  ggplot()+
  geom_point(aes(y = species, x = precision))+
  xlab("gamma estimate")+
  ylab("Method & Insect Group")+
  theme_bw()+
  theme(plot.title = element_text(size = 7),
        axis.text.x = element_text(size = 8), 
        axis.text.y =element_text(size = 5) )

#save figure
ggsave("Figures/interceptLambda.png", 
       plot = ggInterceptLambda, 
       dpi = 320, 
       height = 8, 
       width = 10, 
       units = "cm")

ggsave("Figures/interceptLambdaSD.png", 
       plot = ggInterceptLambdaSD, 
       dpi = 320, 
       height = 8, 
       width = 10, 
       units = "cm")


####################
# Shannon Diversity
#####################

# Getting Shannon Index from data
shannonEstimates <- lapply(allData, extractShan)%>%
  do.call("cbind", .)%>%
  as.data.frame()%>%
  mutate(site = simulations_all[[1]]$sites$X1km_square)

colnames(shannonEstimates) <- c("bbIDM",
                                "bbSPE",
                                "bbGRO", 
                                "hvIDM",  
                                "hvSPE", 
                                "hvGRO", 
                                "sbIDM", 
                                "sbSPE", 
                                "sbGRO", 
                                "sites")

#saving the shannon diversity with the region and site ID
gr_ref <- read_csv("gr_ref.csv")
names(gr_ref)[1] <- "sites"
sites <- data.frame(sites = simulations_all[[1]]$sites$X1km_square)
merged_data_sites <- merge(x=sites, 
                           y=gr_ref, 
                           by="sites", 
                           all.x = TRUE, 
                           all.y = FALSE)%>%
  select(sites, region)%>%
  arrange(sites)

#Fill in the missing region for a location
merged_data_sites[62,2] <- "ENG"
shannonEstimates <- merge(shannonEstimates, 
                            merged_data_sites, 
                            by = "sites", 
                            all.x = TRUE, 
                            all.y = FALSE)

names(shannonEstimates) <- c("sites",
                             "bb_IDM",
                             "bb_SPE",
                             "bb_GRO", 
                             "hv_IDM",  
                             "hv_SPE", 
                             "hv_GRO", 
                             "sb_IDM", 
                             "sb_SPE", 
                             "sb_GRO", 
                             "region")

write.csv(shannonEstimates,"Figures/shannon_index.csv", row.names = FALSE)

# Format the shannon index plot for plot
siteID <- shannonEstimates$sites
shannonEstimatesLongFormat <- pivot_longer(shannonEstimates, 
                             cols = c(-sites, -region),
                             names_to = c("group", "model"),
                             names_sep = "_")

# Mean Shannon indices across all sites
shannonBoxplot <- ggplot(data = shannonEstimatesLongFormat) +
  geom_boxplot(aes(x = group, y = value,
                   fill = model,
                   color = model),
               alpha = 0.2) +
  #facet_wrap(~year) +
  ylab("Mean Shannon index across sites") +
  xlab("Insects") +
  theme_bw()+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                    name="Method", labels = c("Insect group", "IDM", "Species only")) +
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                     name="Method", labels = c("Insect group", "IDM", "Species only")) +
  scale_x_discrete(labels= c("Bumblebees", "Hoverflies", "Solitary bees"))

#save plot
ggsave("Figures/shannonBoxplot.png", 
       plot = shannonBoxplot, 
       dpi = 320, 
       height = 8, 
       width = 10, 
       units = "cm")


# Plot the Shannon Indices on the map
# map of sites (represented at 10Km because too small to see at 1Km)
plot_GIS(UK[UK$REGION == "Great Britain",], xlab="", ylab="", show.axis = FALSE, show.grid = FALSE, 
         fill.col = "lightgrey", no.margin = TRUE)

#plotUK_gr(shannon_long$sites, gr_prec = 10000, border = "black")
#Get coordinates for the sites and add to shannon estimates
coords <- gr_let2num(shannonEstimatesLongFormat$sites)
allDataShannon <- cbind(shannonEstimatesLongFormat, coords)

UK_countries_lowres@data$ID = rownames(UK_countries_lowres@data)
UK_countries_lowres_df <- ggplot2::fortify(UK_countries_lowres)

allDataShannon <- allDataShannon %>%
  mutate(model = case_when(model == "GRO" ~ "Group only",
                           model == "SPE" ~ "Species only",
                           model == "IDM" ~ "IDM"))

shannonMap <- lapply(as.list(c("bb", "hv" , "sb")), function(x){
  ggplot(data = allDataShannon%>%
                           filter(group == x) ) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~ model , nrow = 1, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon", limits = c(0, 4))
})

gg4 <- ggpubr::ggarrange(plotlist = shannonMap,
                         ncol = 1, nrow = 3, 
                         common.legend = TRUE, 
                         legend = "bottom",
                         labels = c("a) Bumblebees",
                                    "b) Hoverflies",
                                    "c) Solitarybees")) 

ggsave("Figures/shannonMap.png", 
       plot = gg4, 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")

############################
# Hills indices
#########################

#Retrieve estimated hills indices
estimatedHillsIndices  <- lapply(allData, function(x){
  lapply(as.list(0:2), function(y){
    hillsIndex(x, y)
  })
})%>%
  purrr::flatten(.)%>%
  do.call("cbind", .)%>%
  as.data.frame()%>%
  mutate(site = simulations_all[[1]]$sites$X1km_square)

colnames(estimatedHillsIndices) <- c("bb_IDM_0", "bb_IDM_1", "bb_IDM_2", 
                                     "bb_SPE_0", "bb_SPE_1", "bb_SPE_2",
                                     "bb_GRO_0", "bb_GRO_1", "bb_GRO_2",
                                     "hv_IDM_0", "hv_IDM_1", "hv_IDM_2",
                                     "hv_SPE_0", "hv_SPE_1", "hv_SPE_2",
                                     "hv_GRO_0", "hv_GRO_1", "hv_GRO_2",
                                     "sb_IDM_0", "sb_IDM_1", "sb_IDM_2",
                                     "sb_SPE_0", "sb_SPE_1", "sb_SPE_2",
                                     "sb_GRO_0", "sb_GRO_1", "sb_GRO_2",
                                     "sites")

estimatedHillsIndicesAll <- merge( estimatedHillsIndices, 
                                   merged_data_sites, 
                                   by = "sites", 
                                   all.x = TRUE, 
                                   all.y = FALSE)

write.csv(estimatedHillsIndicesAll,"Figures/estimatedHillsIndices.csv", row.names = FALSE)

# Format the data into a longer format
hillsLong <- pivot_longer(estimatedHillsIndicesAll, cols = c(-sites, -region),
                           names_to = c("group", "model", "q"),
                           names_sep = "_")

#Boxplot for Hills Indices
hillsBoxplot <- ggplot(data = hillsLong) +
  geom_boxplot(aes(x = group, y = value,
                   fill = model,
                   color = model),
               alpha = 0.2) +
  facet_wrap(~q) +
  ylab("Mean Hills Indices across sites") +
  xlab("Insects") +
  theme_bw()+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                    name="Method", labels = c("Insect group", "IDM", "Species only")) +
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                     name="Method", labels = c("Insect group", "IDM", "Species only")) +
  scale_x_discrete(labels= c("Bumblebees", "Hoverflies", "Solitary bees"))

#save plot
ggsave("Figures/hillsBoxplot.png", 
       plot = hillsBoxplot, 
       dpi = 320, 
       height = 8, 
       width = 10, 
       units = "cm")

#get coordinates for the sites
coords <- gr_let2num(hillsLong$sites)
hillsPlotData <- cbind(hillsLong, coords)

UK_countries_lowres@data$ID = rownames(UK_countries_lowres@data)
UK_countries_lowres_df <- ggplot2::fortify(UK_countries_lowres)

#prepare data for map plot
hillsPlotData <- hillsPlotData %>%
  mutate(model = case_when(model == "GRO" ~ "Group only",
                           model == "SPE" ~ "Species only",
                           model == "IDM" ~ "IDM"))

#Map plot of hoverflies
HillsPlotData <- lapply(as.list(c("bb", "hv" , "sb")), function(x){
  ggplot(data = hillsPlotData %>%
           filter(group == x)) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~q+model, nrow = 3, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Hills Indices")
})

#save the bumblebees
ggsave("Figures/bbHillsPlot.png", 
       plot = HillsPlotData[[1]], 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")

#save the hoverflies
ggsave("Figures/hvHillsPlot.png", 
       plot = HillsPlotData[[2]], 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")

#save the solitarybees
ggsave("Figures/sbHillsPlot.png", 
       plot = HillsPlotData[[3]], 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")
}

DataFormatting("covariate_inter")

##################
# Data Exploration
#############
gr_ref <- read_csv("gr_ref.csv")
names(gr_ref)[1] <- "sites"
sites <- data.frame(sites = simulations_all[[1]]$sites$X1km_square)
merged_data_sites <- merge(x=sites, 
                           y=gr_ref, 
                           by="sites", 
                           all.x = TRUE, 
                           all.y = FALSE)%>%
  select(sites, region)%>%
  arrange(sites)

#load data
load("data_idm.RData")

#Function that retrieves the species occupancy and group count data and plots
dataForRawData <- lapply(simulations_all, function(x){
rawDataX <- apply(x$mat.species, c(1,4), function(y){
  tt <- sum(is.finite(y))
  if(tt != 0){
  ret <- mean(y>0, na.rm = TRUE)
  }else{
    ret <- NA
  }

})

rawDataY <- apply(x$mat.genus, c(1,3), function(y){
  tt <- sum(is.finite(y))
  if(tt != 0){
    ret <- mean(y>0, na.rm = TRUE)
  }else{
    ret <- NA
  }
  
})

colnames(rawDataX) <- c("X2017", "X2018")
colnames(rawDataY) <- c("Y2017", "Y2018")
coords <- gr_let2num(x$sites$X1km_square)

all_data <- cbind(rawDataX, rawDataY, coords)%>%
  melt(., id.vars = c("EASTING", "NORTHING"))
})%>%
  do.call("rbind", .)%>%
  mutate(insectGroup = rep(c("bumblebees", "hoverflies", "solitarybees"), each = 296))



UK_countries_lowres@data$ID = rownames(UK_countries_lowres@data)
UK_countries_lowres_df <- ggplot2::fortify(UK_countries_lowres)

#plot of formatted data
rawDataPlot <- ggplot(data = dataForRawData) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~variable + insectGroup, ncol = 6) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon",  
                       limits = c(0, 1), 
                       breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                                  0.6, 0.7, 0.8, 0.9, 1),
                       na.value = "red")
ggsave("Figures/rawDataPlot.png", 
       plot = rawDataPlot, 
       dpi = 320, 
       height = 10, 
       width = 20, 
       units = "cm")
