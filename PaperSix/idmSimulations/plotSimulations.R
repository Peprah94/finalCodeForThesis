library(readr)
library(dplyr)
library(ggplot2)
library(svglite)

shanEstimates0 <- read_csv("simEstimates/shanSimulationsEstimates0.csv")
shanEstimates1 <- read_csv("simEstimates/shanSimulationsEstimates1.csv")
shanEstimates2 <- read_csv("simEstimates/shanSimulationsEstimates2.csv")
shanEstimates3 <- read_csv("simEstimates/shanSimulationsEstimates3.csv")
shanEstimates4 <- read_csv("simEstimates/shanSimulationsEstimates4.csv")
shanEstimates5 <- read_csv("simEstimates/shanSimulationsEstimates5.csv")

shanEstimates <- rbind(shanEstimates0,
                       shanEstimates1,
                       shanEstimates2,
                       shanEstimates3,
                       shanEstimates4,
                       shanEstimates5)

shanEst <- shanEstimates%>%
  ggplot()+
  geom_boxplot(aes(x = model, y = V1, fill = truth))+
 # facet_grid(~truth)+
  ylab("MAE(H')")+
  theme_article()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 25))

#text = element_text(size = 20),

# save shan est
ggsave(filename = "Figures/simHills.png",
       shanEst,
       width = 10,
       height = 6,
       dpi = 320)

parEstimates0 <- read_csv("simEstimates/parsSimulationsEstimates0.csv")
parEstimates1 <- read_csv("simEstimates/parsSimulationsEstimates1.csv")
parEstimates2 <- read_csv("simEstimates/parsSimulationsEstimates2.csv")
parEstimates3 <- read_csv("simEstimates/parsSimulationsEstimates3.csv")
parEstimates4 <- read_csv("simEstimates/parsSimulationsEstimates4.csv")
parEstimates5 <- read_csv("simEstimates/parsSimulationsEstimates5.csv")

parEstimates <- rbind(parEstimates0,
                       parEstimates1,
                      parEstimates2,
                      parEstimates3,
                      parEstimates4,
                      parEstimates5)%>%
  as.data.frame()%>%
  mutate(muSpecies = if_else(model %in% c("GCMCO", "GCMSH"), NA, muSpecies),
         sigmaSpecies = if_else(model %in% c("GCMCO", "GCMSH"), NA, sigmaSpecies),
         sigmaAlphaSpecies = if_else(model %in% c("GCMCO", "GCMSH"), NA, sigmaAlphaSpecies),
         sigmaAlphaSites = if_else(model %in% c("GCMCO", "GCMSH"), NA, sigmaAlphaSites),
         sigmaLambda = if_else(model %in% c("GCMSH", "IDMSH", "SOM"), NA, sigmaLambda))


colnames(parEstimates)

parsEst <- parEstimates%>%
  mutate(muLatitude = muLatitude + 0.5,
         muSpecies = muSpecies +3,
         sigmaSpecies = sigmaSpecies -1,
         sigmaSites = sigmaSites -1,
         #sigmaLatitude = sigmaLatitude -1,
         sigmaLambda = sigmaLambda - 1,
         sigmaAlphaSites = sigmaAlphaSites - 1,
         sigmaAlphaSpecies = sigmaAlphaSpecies - 1)%>%
  melt(id.vars = c("model", "truth"))%>%
  ggplot()+
  geom_boxplot(mapping = aes(x = variable, y = value, fill = model))+
  ylim(c(-1, 1.5))+
  facet_grid(~truth) + ylab("Bias of Parameters")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme_article()+
  theme(axis.text  = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom")

parsEst

ggsave(filename = "Figures/simBiasParsEsts.png",
       parsEst,
       width = 16,
       height = 8,
       dpi = 320)
