library(readr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
library(grid)
library(egg)
library(ggh4x)
theme_set(theme_classic(base_size = 14))

allData <- read_csv("formatDataResult/allDataPutTogether.csv")

#text_high <- textGrob("Highest\nvalue", gp=gpar(fontsize=13, fontface="bold"))
#text_low <- textGrob("Lowest\nvalue", gp=gpar(fontsize=13, fontface="bold"))



#new facet labels for estNval
estNValLabs <- c("(x) 10", "(y) 200", "(z) 500")
names(estNValLabs) <- c("10", "200", "500")
trueMisclassIncreaseLabs <- c("(a) Baseline", "(b) Decrease")
names(trueMisclassIncreaseLabs) <- c("0", "6")
simMethodLabs <- c("(ii) full", "(iii) reduced", "(i) correlation")
names(simMethodLabs) <- c("full", "reduced", "correlation")
fills <- c("#66CC00","#66CC00","#66CC00","#66FF33","#66FF33","#66FF33","#00FFFF","#00CCCC","#009999")


predPerformance1 <- allData%>%
  select(3, 5, 6, 12, 15, 8, 13)%>%
  melt(., id.vars = c("simMethod", "estMethod", "estNval", "trueMisclassIncrease"))%>%
  filter(estNval %in% c("200"))%>%
  filter(!trueMisclassIncrease %in% c("4", "2"))%>%
  filter(!variable %in% c("precision"))%>%
  mutate(estMethod = ifelse(estMethod == "onlyCov", "fixed-covariate",
                            ifelse(estMethod == "fixedInterCov", "fixed-intercov", estMethod)))%>%
  ggplot( )+
  geom_boxplot(mapping = aes(x = as.factor(estMethod), y = value, fill = as.factor(variable)), outlier.alpha = 0)+
  xlab("")+
  ylab("Accuracy and recall")+
  #labs(title = "Change in number of misclassified samples in validation samples")+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_fill_manual(values = c("#E69F00", "#F0E442"))+ #"#D55E00"
  facet_grid( trueMisclassIncrease  ~  simMethod, scales = "free_y",
              labeller = labeller(#estNval = estNValLabs,
                                  trueMisclassIncrease = trueMisclassIncreaseLabs,
                                  simMethod = simMethodLabs))+
  #scale_y_continuous(sec.axis = sec_axis(~.+10, name='Number of validation samples')) +
  coord_trans()+
 # theme_classic()+
  theme(legend.position = "right", legend.title = element_blank(),
        #strip.placement = 'outside',
        strip.background = element_blank(),
       # plot.title = element_text(hjust = 0.5, size = 12, colour = "red"),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        text = element_text(size = 20),
       axis.text.x = element_blank()
        #axis.title.y.right = element_text(hjust = 0.5, size = 12, colour = "red")
       )

predPerformance1

predPerformance2 <- allData%>%
  select(3, 5, 6, 12, 15, 8, 13)%>%
  melt(., id.vars = c("simMethod", "estMethod", "estNval", "trueMisclassIncrease"))%>%
  filter(estNval %in% c("200"))%>%
  filter(!trueMisclassIncrease %in% c("4", "2"))%>%
  filter(variable %in% c("precision"))%>%
  mutate(estMethod = ifelse(estMethod == "onlyCov", "fixed-covariate",
                            ifelse(estMethod == "fixedInterCov", "fixed-intercov", estMethod)))%>%
  ggplot( )+
  geom_boxplot(mapping = aes(x = as.factor(estMethod), y = value, fill = as.factor(variable)), outlier.alpha = 0)+
  xlab("Model")+
  ylab("Precision")+
  #labs(title = "Change in number of misclassified samples in validation samples")+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_fill_manual(values = c("#D55E00", "#F0E442"))+ #"#D55E00"
  facet_grid( trueMisclassIncrease ~  simMethod, scales = "free_y",
              labeller = labeller(#estNval = estNValLabs,
                trueMisclassIncrease = trueMisclassIncreaseLabs,
                simMethod = simMethodLabs))+
  #scale_y_continuous(sec.axis = sec_axis(~.+10, name='Number of validation samples')) +
  coord_trans()+
  # theme_classic()+
  theme(legend.position = "right", legend.title = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(),
        #strip.placement = 'outside',
        #strip.background = element_rect(colour = "black", fill = "white"),
        #plot.title = element_text(hjust = 0.5, size = 12, colour = "red"),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.text.y.right = element_blank(),
       text = element_text(size = 20),
        axis.title.y.right = element_text(hjust = 0.5, size = 12, colour = "red"))

predPerformance2

performance <- ggpubr::ggarrange(predPerformance1, predPerformance2,
                  ncol = 1,
                  common.legend = TRUE,
                  labels = c("A)" , "B)"))

ggsave(filename = "metrics.png", plot = performance, dpi = 700, width = 30, height = 30, units = "cm")

#Summarise bias in parameters for table
#Bias in parameters
estimateBiasPars <- allData%>%
 dplyr::select("simMethod", "estMethod", "p","waic" ,"trueMisclassIncrease", "estNval",
               "biasBeta01"  ,         "biasBeta11"  ,         "biasBeta12" ,
               "biasGamma011"  ,       "biasGamma021"  ,       "biasGamma012" ,
                "biasGamma022" ,        "biasGamma111"  ,       "biasGamma121",
                "biasGamma112" ,        "biasGamma122")%>%
  filter(estNval %in% c("200"))%>%
  filter(!trueMisclassIncrease %in% c("4", "2"))%>%
  group_by(simMethod, estMethod, trueMisclassIncrease) %>%
  mutate(estMethod = ifelse(estMethod == "onlyCov", "fixed-covariate",
                            ifelse(estMethod == "fixedInterCov", "fixed-intercov", estMethod)))%>%
  summarise(p = paste0(round(mean(p),4),"(", round(sd(p),4), ")"),
            waic = paste0(round(mean(waic),4),"(", round(sd(waic),4), ")"),
            #sdP = sd(p),
            biasBeta01 = paste0(round(mean(biasBeta01),4),"(", round(sd(biasBeta01),4), ")"),
            biasBeta11 = paste0(round(mean(biasBeta11),4),"(", round(sd(biasBeta11),4), ")") ,
            biasBeta12 = paste0(round(mean(biasBeta12),4),"(", round(sd(biasBeta12),4), ")"),
          biasGamma011= paste0(round(mean(biasGamma011),4),"(", round(sd(biasGamma011),4), ")")  ,
          biasGamma021 = paste0(round(mean(biasGamma021),4),"(", round(sd(biasGamma021),4), ")"),
          biasGamma012 = paste0(round(mean(biasGamma012),4),"(", round(sd(biasGamma012),4), ")"),
            biasGamma022 = paste0(round(mean(biasGamma022),4),"(", round(sd(biasGamma022),4), ")"),
          biasGamma111 = paste0(round(mean(biasGamma111),4),"(", round(sd(biasGamma111),4), ")") ,
          biasGamma121 = paste0(round(mean(biasGamma121),4),"(", round(sd(biasGamma121),4), ")"),
            biasGamma112 = paste0(round(mean(biasGamma112),4),"(", round(sd(biasGamma112),4), ")"),
          biasGamma122 = paste0(round(mean(biasGamma122),4),"(", round(sd(biasGamma122),4), ")"))%>%
  ungroup()%>%
  mutate(trueMisclassIncrease = ifelse(trueMisclassIncrease == 0, "Baseline", "Decreased"))

write.csv(estimateBiasPars, file = paste0("formatDataResult/estimateBiasPars.csv"))


