#load data
source("idmFunctions/functionForEstimation.R")
# cross Validation

cvRes <- lapply(as.list(c("bumblebees", "hoverflies", "solitarybees")), function(x){
load(paste0("crossValidation/",x,"/estimateCrossValidate.RData"))
idm <- estimatesCaseStudy
})

load("/Volumes/kwakupa/idmProject/Figures/correlationMatrices.RData")

title_text <- c("a) IDMSH",
                "b) SOM",
                "c) IDMCO",
                #"b) bumblebees \n Group only",
               # "c) bumblebees \n Species only",
               # "c) bumblebees \n Species only",


                "a) IDMSH",
                "b) SOM",
                "c) IDMCO",
                #"e) hoverflies \n Group only",
               # "f) hoverflies \n Species only",
                #"f) hoverflies \n Species only",


                "a) IDMSH",
                "b) SOM",
                "c) IDMCO"#,
                #"h) solitarybees \n Group only",
                #"i) solitarybees \n species only",
                #"i) solitarybees \n species only"
               )

#Plot the correlation matrices
ggCorrplots <-  lapply(seq_along(correlationMatrices), function(x){
  ggcorrplot(correlationMatrices[[x]],
             hc.order = TRUE,
             type = "upper",
             tl.cex = 5,
             tl.srt = 90,
             lab_size = 5,
             title = title_text[[x]])+
    theme_classic() +
    theme(plot.title = element_text(size = 9),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 5) )+
    ylab(" ")+
    xlab(" ")

})

# Arrange the plots
#$shared model
gg1 <- ggpubr::ggarrange(ggCorrplots[[1]], ggCorrplots[[2]], ggCorrplots[[3]],
                         ncol = 2, nrow = 2, common.legend = TRUE,
                         legend = "right")

gg2 <- ggpubr::ggarrange(ggCorrplots[[4]], ggCorrplots[[5]], ggCorrplots[[6]],
                         ncol = 2, nrow = 2, common.legend = TRUE,
                         legend = "right")

gg3 <- ggpubr::ggarrange(ggCorrplots[[7]], ggCorrplots[[8]], ggCorrplots[[9]],
                         ncol = 2, nrow = 2, common.legend = TRUE,
                         legend = "bottom")

#covariate model
# gg3 <- ggpubr::ggarrange(ggCorrplots[[2]], ggCorrplots[[4]], ggCorrplots[[6]],
#                          ggCorrplots[[8]], ggCorrplots[[10]], ggCorrplots[[12]],
#                          ggCorrplots[[14]], ggCorrplots[[16]], ggCorrplots[[18]],
#                          ncol = 3, nrow = 3, common.legend = TRUE,
#                          legend = "right",
#                          heights = c(1.5,2,2))

figure <- ggpubr::annotate_figure(gg1, left = "Species list", bottom = "Species list")
figure1 <- ggpubr::annotate_figure(gg2, left = "Species list", bottom = "Species list")
figure2 <- ggpubr::annotate_figure(gg3, left = "Species list", bottom = "Species list")
#figure1 <- ggpubr::annotate_figure(gg3, left = "Species list", bottom = "Species list", top = "Covariate model")

#save plot
ggsave("Figures/speciesAssociationBumblebees.png", plot = figure, dpi = 320,
       height = 18, width = 18, units = "cm")

ggsave("Figures/speciesAssociationHoverflies.png", plot = figure1, dpi = 320,
       height = 18, width = 18, units = "cm")

ggsave("Figures/speciesAssociationSolitarybees.png", plot = figure2, dpi = 320,
       height = 18, width = 18, units = "cm")
#ggsave("Figures/speciesAssociationCovariate.png", plot = figure1, dpi = 320,
#       height = 18, width = 18, units = "cm")

