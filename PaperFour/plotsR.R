#plot of example 1
library(dplyr)
library(mcmcse)
library(ggplot2)
library(ggpubr)
library(INLA)
library(ggmcmc)

fill.colors <- c("INLA" = "darkgreen",
                 "iNim-RW" = "#FF3300",
                 #"MCMC" = "#00CCFF",
                 #"iNIM2-RW" = "#E69F00",
                 "MCMC" = "#0033FF",
                 "iNim2-RW" = "#E69F00"
                 #"iNim2-RW" = "#FF6600"#,
                 #"truth" = "black"
                 )
#########################
# Bivariate regression
#########################

load("bivariateRegression/inlaWithNimbleAlt/inlaNimBivariateRegressionAlt.RData")
model4 <- inlaNimBivariateAlt[[1]]
load("bivariateRegression/inlaWithNimble/inlaNimBivariateRegression1.RData")
model1 <- inlaNimBivariate
load("bivariateRegression/inlaWithNimble/inlaNimBivariateRegression2.RData")
model2 <- inlaNimBivariate
load("bivariateRegression/inlaWithNimble/inlaNimBivariateRegression3.RData")
model3 <- inlaNimBivariate
#fit INLA model
load("bivariateRegression/data_for_simulations.RData")

# set the data and the offset
data <- data.frame(y = bivariateSims$y,
                   x1 = bivariateSims$x[,1],
                   x2 = bivariateSims$x[,2])

res = INLA::inla(y ~ 1 + x1 + x2,
                 data = data,
                 #family = family,
                 verbose=FALSE,
                 control.fixed = list(prec.intercept = 0.001),
                 control.predictor = list(compute = TRUE),
                 control.compute=list(config = TRUE))

# conditional marginal log-likelihood
samples <- inla.posterior.sample(5000, res)
sampleshyper <- inla.hyperpar.sample(5000, res)
inlaSamples <- lapply(samples, function(x){
  c(c(tail(x[[2]], 3)))#, c(sampleshyper), "INLA")
})%>%
  do.call("rbind",.)%>%
  as.data.frame()%>%
  mutate(tau = c(sampleshyper),
         method = "INLA")

colnames(inlaSamples) <- c("a", "beta[1]",  "beta[2]",  "tau"  ,   "method" )


#extract samples
method <- c("iNim-RW", "iNim-AFSS" , "MCMC")
ret <- rbind(cbind(model1$inlamcmc$mcmc.out$samples, method[1]),
             cbind(model2$inlamcmc$mcmc.out$samples, method[2]),
             cbind(model3$mcmc$mcmc.out$samples, method[3]))%>%
  as.data.frame()

colnames(ret)[5] <- "method"
names <- colnames(ret)

colnames(ret) <- c("a","beta[1]",  "beta[2]",  "tau"  ,   "method" )

#extract samples
method <- c("iNim2-RW" )
ret1 <- cbind(model4$mcmc.out$samples, method)
colnames(ret1) <- c("beta[1]",  "beta[2]",  "a","tau"  ,   "method" )

columns_to_plot <- c("beta[1]",  "beta[2]",  "a","tau")

allResults <- rbind(ret, ret1, inlaSamples)%>%
  as.data.frame()%>%
  mutate(method = as.factor(method))%>%
  mutate_at(vars(all_of(columns_to_plot)), as.numeric)

plotList <- list()
true_vals <- c(3, -3, 2, 1)
for(i in seq_along(columns_to_plot)){

  col_greek <- parse(text = paste(columns_to_plot[i]))

  plot <- allResults%>%
    filter(., !method %in% c("iNim-AFSS"))%>%
    ggplot(., aes(x = .data[[columns_to_plot[i]]], linetype = as.factor(method), col = as.factor(method)))+
  geom_density()+
    labs(title = col_greek)+
    xlab("")+
    theme_bw()+
    geom_vline(xintercept = true_vals[i])+
    scale_color_manual(values = fill.colors)+
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.5, "cm"))
    #ggpubr::theme_classic2()
  #print(plot)

  plotList[[i]] <- plot
}

contPlot <- allResults%>%
  dplyr::filter(!method %in% c("INLA", "iNim-AFSS"))%>%
  ggplot(., mapping = aes(x = `beta[1]`,
                          y = `beta[2]`))+
  geom_density_2d()+
  geom_point(mapping = aes(x = 3, y = -3), shape = 4, size = 3)+
  xlab(parse(text = paste(columns_to_plot[1])))+
  ylab(parse(text = paste(columns_to_plot[2])))+
  facet_wrap(~method, nrow = 2, ncol = 2, scales = "free")+
  theme_bw()

ggsave("results/bivariateContour.png",
       plot = contPlot,
       width = 7,
       height = 4,
       units = "in")

bivPlot <- ggpubr::ggarrange(plotlist = plotList,
                  ncol = 2,
                  nrow = 2,
                  common.legend = TRUE,
                  legend = "bottom")
ggsave("results/bivariateRegression.png",
       plot = bivPlot,
       width = 7,
       height = 4,
       units = "in")

## Plot for tau
retTauInlaRW <- ggs(model1$inlamcmc$mcmc.out$samples)%>%
  ggs_traceplot(., family = "sigma")+
  ggtitle("iNim-RW")

retTauInlaMCMC <- ggs(model3$mcmc$mcmc.out$samples)%>%
  ggs_traceplot(., family = "sigma")+
  ggtitle("MCMC")

retTauInlaRW2 <- ggs(model4$mcmc.out$samples)%>%
  ggs_traceplot(., family = "sigma")+
  ggtitle("MCMC")

traceplots <- ggarrange(retTauInlaRW,
                        retTauInlaMCMC,
                        retTauInlaRW2,
                        ncol = 1)
ggsave("results/bivariateRegressionSigma.png",
       plot = traceplots,
       width = 7,
       height = 12,
       units = "in")             

#Estimate Effective sample sizes and time taken
method <- c("iNim-RW", "iNim-AFSS" , "MCMC")
ret <- rbind(c(ess(model1$inlamcmc$mcmc.out$samples), method[1], as.numeric(model1$inlamcmc$timeTaken)/3600, c(ess(model1$inlamcmc$mcmc.out$samples)/as.numeric(model1$inlamcmc$timeTaken))),
             c(ess(model2$inlamcmc$mcmc.out$samples), method[2], as.numeric(model2$inlamcmc$timeTaken)/3600, c(ess(model2$inlamcmc$mcmc.out$samples)/as.numeric(model2$inlamcmc$timeTaken))),
             c(ess(model3$mcmc$mcmc.out$samples), method[3], model3$mcmc$timeTaken, c(ess(model3$mcmc$mcmc.out$samples)/as.numeric(model3$mcmc$timeTaken)))
             )%>%
  as.data.frame()

colnames(ret)[5] <- "method"
colnames(ret)[6] <- "timeTaken"
names <- colnames(ret)
colnames(ret) <- c("a_ESS","beta[1]_ESS",  "beta[2]_ESS",  "tau_ESS"  ,   "method", "timeTaken","a_Eff","beta[1]_Eff",  "beta[2]_Eff",  "tau_Eff"  )

method <- c("iNim2-RW" )
rr <- c(ess(model4$mcmc.out$samples))
rrEst <- c(ess(model4$mcmc.out$samples))/as.numeric(model4$timeTaken)
ret1 <- data.frame(rr[3], rr[1], rr[2], rr[4], method, as.numeric(model4$timeTaken)/3600, rrEst[3], rrEst[1], rrEst[2], rrEst[4])
colnames(ret1) <- colnames(ret) <- c("a_ESS","beta[1]_ESS",  "beta[2]_ESS",  "tau_ESS"  ,   "method", "timeTaken","a_Eff","beta[1]_Eff",  "beta[2]_Eff",  "tau_Eff"  )

res <- rbind(ret,
             ret1)%>%
  filter(!method %in% c("iNim-AFSS"))
# save results
write.csv(res, file = "results/bivariateESS.csv")

#### Bayesian lasso regression

#load results
load("bayesianLasso/inlaWithNimble/inlaNimBayesianLasso1.RData")
model1 <- bayesianLasso
load("bayesianLasso/inlaWithNimble/inlaNimBayesianLasso2.RData")
model2 <- bayesianLasso
load("bayesianLasso/inlaWithNimble/inlaNimBayesianLasso3.RData")
model3 <- bayesianLasso
load("bayesianLasso/inlaWithNimbleAlt/inlaNimBayesianLassoAlt.RData")
model4 <- bayesianLassoAlt

#get results from lasso regression
load("bayesianLasso/hitters_data.RData")
library(glmnet)
library(ISLR)
df <- lassoDataDF
x <- df$x
y <-df$y
n.beta <- ncol(df$x)

# ml estimates
ml = summary(lm(y~-1 + x, data = df))$coefficients[,1:2]

#Indices for train/test model
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)

#Grid for lambda parameter in lasso
grid <- 10^seq(10, -2, length = 100)

#Fit lasso model for several values of lambda
lasso.mod <- glmnet(x[train, ] , y[train], alpha = 1, lambda = grid,intercept = F)

#CV
set.seed(1)
cv.out <- cv.glmnet(x[train, ], y[train], alpha = 1,intercept=F)

#Take best lambda for lasso model
bestlam <- cv.out$lambda.min

#Predcit with lasso on test data
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test, ])

#Fit model to complete dataset
out <- glmnet(x, y, alpha = 1, lambda = grid,intercept=F)
lasso.coef <- predict(out, type = "coefficients", s = bestlam)

method <- c("iNim-RW", "iNim-AFSS" , "MCMC")
ret <- rbind(cbind(model1$inlamcmc$mcmc.out$samples, method[1]),
             cbind(model2$inlamcmc$mcmc.out$samples, method[2]),
             cbind(model3$mcmc$mcmc.out$samples, method[3]))%>%
  as.data.frame()

colnames(ret)[8] <- "method"
names <- colnames(ret)

#extract samples
method <- c("iNim2-RW" )
ret1 <- cbind(model4[[1]]$mcmc.out$samples, method)
#colnames(ret1) <- c("beta[1]",  "beta[2]",  "a","sigma"  ,   "method" )

columns_to_plot <- c("beta[1]",  "beta[2]", "beta[3]", "beta[4]", "beta[5]" )

allResults <- rbind(ret, ret1)%>%
  as.data.frame()%>%
  mutate(method = as.factor(method))%>%
  mutate_at(vars(all_of(columns_to_plot)), as.numeric)

plotList <- list()

names_plot <- c("AtBat",  "Hits", "HmRun", "Runs", "RBI" )
lassoVals <- as.numeric(lasso.coef)[-1]
for(i in seq_along(columns_to_plot)){

  col_greek <- parse(text = paste(columns_to_plot[i]))

  plot <- allResults%>%
    dplyr::filter(!method %in% c("iNim-AFSS"))%>%
    ggplot(., aes(x = .data[[columns_to_plot[i]]], linetype = as.factor(method), col = as.factor(method)))+
    geom_density()+
    labs(title = names_plot[i])+
    xlab("")+
    theme_bw()+
    geom_vline(xintercept = lassoVals[i])+
    scale_color_manual(values = fill.colors)+
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.5, "cm"))
  #print(plot)

  plotList[[i]] <- plot
}

lassoPlot <- ggpubr::ggarrange(plotlist = plotList,
                  ncol = 2,
                  nrow = 3,
                  common.legend = TRUE,
                  legend = "bottom")


ggsave("results/lassoRegression.png",
       plot = lassoPlot,
       width = 7,
       height = 6,
       units = "in")



# Extract estimates for table
lasso.fitted <- predict(out, s = bestlam, newx = x)

method <- c("iNim-RW", "iNim-AFSS" , "MCMC")
ret <- rbind(c(model1$inlamcmc$mcmc.out$summary[,1], method[1]),
             c(model2$inlamcmc$mcmc.out$summary[,1], method[2]),
             c(model3$mcmc$mcmc.out$summary[,1], method[3]))%>%
  as.data.frame()

retSD <- rbind(c(model1$inlamcmc$mcmc.out$summary[,3], method[1]),
               c(model2$inlamcmc$mcmc.out$summary[,3], method[2]),
               c(model3$mcmc$mcmc.out$summary[,3], method[3]))%>%
  as.data.frame()

colnames(ret)[8] <- "method"
names <- colnames(ret)

#extract samples
method <- c("iNim2-RW")
ret1 <- c(bayesianLassoAlt[[1]]$mcmc.out$summary[,1], method)
#ret1 <- c(missingCovsAlt[[1]]$mcmc.out$summary[,1], method)
retSD1 <- c(bayesianLassoAlt[[1]]$mcmc.out$summary[,3], method)

resMean <- cbind(rbind(ret,
                       ret1),
                 "mean")
colnames(resMean)[9] <- "metrics"

resSD <- cbind(rbind(retSD,
                     retSD1),
               "sd")

colnames(resSD) <- colnames(resMean)

allResults <- rbind(resMean,
                    resSD)%>%
  filter(!method %in% c("iNim-AFSS"))

# save results
write.csv(allResults, file = "results/bayesianLasso.csv")

##################
# Missing Covariates
#################
load("missingCovariates/inlaWithNimble/inlaNimMissingCovariates1.RData")
model1 <- missingCovs
load("missingCovariates/inlaWithNimble/inlaNimMissingCovariates2.RData")
model2 <- missingCovs
load("missingCovariates/inlaWithNimble/inlaNimMissingCovariates3.RData")
model3 <- missingCovs
load("missingCovariates/inlaWithNimbleAlt/inlaNimMissingCovariatesAlt.RData")
model4 <- missingCovsAlt

#extract samples
method <- c("iNim-RW", "iNim-AFSS" , "MCMC")
ret <- rbind(cbind(model1$inlamcmc$mcmc.out$samples, method[1]),
             cbind(model2$inlamcmc$mcmc.out$samples, method[2]),
             cbind(model3$mcmc$mcmc.out$samples, method[3]))%>%
  as.data.frame()

colnames(ret)[15] <- "method"
names <- colnames(ret)

#extract samples
method <- c("iNim2-RW")
ret1 <- cbind(missingCovsAlt$mcmc.out$samples, method)
colnames(ret1) <- colnames(ret)

columns_to_plot <- c("eta[1]",  "eta[2]", "eta[3]","eta[4]",
                     "eta[5]",  "eta[6]", "eta[7]","eta[8]",
                     "eta[9]")
names_plots <- c("Observation 1",
                 "Observation 3",
                 "Observation 4",
                 "Observation 6",
                 "Observation 10",
                 "Observation 11",
                 "Observation 12",
                 "Observation 16",
                 "Observation 21")

allResults <- rbind(ret, ret1)%>%
  as.data.frame()%>%
  mutate(method = as.factor(method))%>%
  mutate_at(vars(all_of(columns_to_plot)), as.numeric)

plotList <- list()

for(i in seq_along(columns_to_plot)){

  col_greek <- parse(text = paste(columns_to_plot[i]))

  plot <- allResults%>%
    filter(!method %in% c("iNim-AFSS"))%>%
  ggplot(., aes(x = .data[[columns_to_plot[i]]], linetype = as.factor(method), col = as.factor(method)))+
    geom_density()+
    labs(title = names_plots[i])+
    xlab("")+
    theme_bw()+
    scale_color_manual(values = fill.colors)+
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.5, "cm"))
  #print(plot)

  plotList[[i]] <- plot
}


missingPlot <- ggpubr::ggarrange(plotlist = plotList,
                  ncol = 3,
                  nrow = 3,
                  common.legend = TRUE,
                  legend = "bottom")

ggsave("results/missingRegression.png",
       plot = missingPlot,
       width = 9,
       height = 6,
       units = "in")

method <- c("iNim-RW", "iNim-AFSS" , "MCMC")
ret <- rbind(c(model1$inlamcmc$mcmc.out$summary[,1], method[1]),
             c(model2$inlamcmc$mcmc.out$summary[,1], method[2]),
             c(model3$mcmc$mcmc.out$summary[,1], method[3]))%>%
  as.data.frame()

retSD <- rbind(c(model1$inlamcmc$mcmc.out$summary[,3], method[1]),
               c(model2$inlamcmc$mcmc.out$summary[,3], method[2]),
               c(model3$mcmc$mcmc.out$summary[,3], method[3]))%>%
  as.data.frame()

colnames(ret)[15] <- "method"
names <- colnames(ret)

#extract samples
method <- c("iNim2-RW")
ret1 <- c(missingCovsAlt$mcmc.out$summary[,1], method)
retSD1 <- c(missingCovsAlt$mcmc.out$summary[,3], method)
#colnames(ret1) <- colnames(retSD1) <- "method"
#save results
resMean <- cbind(rbind(ret,
                 ret1),
                 "mean")
colnames(resMean)[16] <- "metrics"

resSD <- cbind(rbind(retSD,
               retSD1),
               "sd")

colnames(resSD) <- colnames(resMean)

allResults <- rbind(resMean,
                    resSD)%>%
  filter(!method %in% c("iNim-AFSS"))

# save results
write.csv(allResults, file = "results/missingCovs.csv")


#########################
#   N-Mixture model
###########

library(unmarked)
data("mallard")
library(myphdthesis)
library(nimble)
library(INLA)
library(inlabru)


length <- mallard.site[ , "length"]
elev <- mallard.site[, "elev"]
forest <- mallard.site[, "forest"]
mean.ivel <- mallard.obs$ivel
mean.date <- mallard.obs$date

# unmarked
mallard.umf <- unmarkedFramePCount(y = mallard.y,
                                   siteCovs = mallard.site,
                                   obsCovs = mallard.obs)

out.unmk.2 <- pcount(~ ivel + date + I(date^2) ~ length + elev + forest,
                     mixture = "P",
                     data = mallard.umf)
summary(out.unmk.2)

#extract coefficnents
unmarkedCoef <- coef(out.unmk.2)
unmarkedSE <- SE(out.unmk.2)
names(unmarkedCoef) <- names(unmarkedSE) <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]",
                                              "alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]")

#estimate total abundance with unmarked package
est <- bup(ranef(out.unmk.2))
unmarkedNest <- round(sum(est))
#unmarkedNest <- round(sum(est))

allUnmarkedResults <-  matrix(c(unmarkedCoef, unmarkedNest, "MLE"), nrow =1)%>%
  as.data.frame()

allUnmarkedResultsSE <-  matrix(c(unmarkedSE, NA, "MLE"), nrow =1)%>%
  as.data.frame()

colnames(allUnmarkedResults) <- colnames(allUnmarkedResultsSE) <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]",
                                  "alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]", "Ntotal", "method")

sortedNames <- sort(c("beta[1]", "beta[2]", "beta[3]", "beta[4]",
                      "alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]", "Ntotal", "method"))

allUnmarkedResults <- allUnmarkedResults%>%
  select(all_of(sortedNames))

allUnmarkedResultsSE <- allUnmarkedResultsSE%>%
  select(all_of(sortedNames))

# Load data
load("binomialNMix/inlaWithNimble/binomialNMix1.RData")
model1 <- bayesianNMix

load("binomialNMix/inlaWithNimble/binomialNMix2.RData")
model2 <- bayesianNMix

load("binomialNMix/inlaWithNimble/binomialNMix3.RData")
model3 <- bayesianNMix





method <- c("iNim-RW", "MCMC")
ret <- rbind(c(model1$inlamcmc$mcmc.out$summary[,2], method[1]),
             #c(model2$inlamcmc$mcmc.out$summary[,1], method[2]),
             c(model3$mcmc$mcmc.out$summary[,2], method[2]))%>%
  as.data.frame()%>%
  select(c(192:201))

retSD <- rbind(c(model1$inlamcmc$mcmc.out$summary[,3], method[1]),
               #c(model2$inlamcmc$mcmc.out$summary[,3], method[2]),
               c(model3$mcmc$mcmc.out$summary[,3], method[2]))%>%
  as.data.frame()%>%
  select(c(192:201))

colnames(ret)[10] <- colnames(retSD)[10]  <- "method"
names <- colnames(ret)

ret <- ret%>%
  select(all_of(sortedNames))

retSD <- retSD%>%
  select(all_of(sortedNames))

resMean <- cbind(rbind(ret,
                       allUnmarkedResults),
                 "mean")
colnames(resMean)[11] <- "metrics"

resSD <- cbind(rbind(retSD,
                     allUnmarkedResultsSE),
               "sd")
colnames(resSD)[11] <- "metrics"
#colnames(resSD) <- colnames(resMean)

allResults <- rbind(resMean,
                    resSD)

# save results
write.csv(allResults, file = "results/binomialNmix.csv")

#extract samples
method <- c("iNim-RW",  "MCMC")
ret <- rbind(cbind(model1$inlamcmc$mcmc.out$samples, method[1]),
             #cbind(model2$inlamcmc$mcmc.out$samples, method[2]),
             cbind(model3$mcmc$mcmc.out$samples, method[2]))%>%
  as.data.frame()


colnames(ret)[201] <- "method"
names <- colnames(ret)
#Check marginal distributions
columns_to_plot <- c("beta[1]",  "beta[2]", "beta[3]", "beta[4]", "alpha[1]",
                     "alpha[2]", "alpha[3]", "alpha[4]")

allResults <- rbind(ret)%>%
  as.data.frame()%>%
  mutate(method = as.factor(method))%>%
  mutate_at(vars(all_of(columns_to_plot)), as.numeric)

plotList <- list()

#names_plot <- c("AtBat",  "Hits", "HmRun", "Runs", "RBI" )
#lassoVals <- as.numeric(lasso.coef)[-1]

for(i in seq_along(columns_to_plot)){

  col_greek <- parse(text = paste(columns_to_plot[i]))

  plot <- ggplot(allResults, aes(x = .data[[columns_to_plot[i]]], linetype = as.factor(method), col = as.factor(method)))+
    geom_density()+
    labs(title = col_greek)+
    xlab("")+
    theme_bw()+
    geom_vline(xintercept = unmarkedCoef[i])+
    scale_color_manual(values = fill.colors)+
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.5, "cm"))
  #print(plot)

  plotList[[i]] <- plot
}

binNmixPlot <- ggpubr::ggarrange(plotlist = plotList,
                               ncol = 3,
                               nrow = 3,
                               common.legend = TRUE,
                               legend = "bottom")


ggsave("results/binNmix.png",
       plot = binNmixPlot,
       width = 7,
       height = 6,
       units = "in")


library(ggmcmc)

ggs(model1$inlamcmc$mcmc.out$samples)%>%
  ggs_traceplot(., family = "alpha")

##############
# Zero-inflated Poisson
##############
library(pscl)
load("zip/zinb.RData")
summary(zinb)

#Maximum likehood estimates
# Full data
d <- zinb

# ZIP
ml.res <- summary(zeroinfl(count ~ child + camper | persons, data = d))
coefs <- coef(ml.res)
coefsEst <- c(coefs$count[,1], coefs$zero[,1], "MLE")
coefsSD <- c(coefs$count[,2], coefs$zero[,2], "MLE")


load("/Volumes/kwakupa/inlamcmc/zip/inlaWithNimble/inlaNimZIPSmall1.RData")
model1 <- zipModel

load("/Volumes/kwakupa/inlamcmc/zip/inlaWithNimble/inlaNimZIPSmall3.RData")
model3 <- zipModel

load("/Volumes/kwakupa/inlamcmc/zip/inlaWithNimbleAlt/zipModelAltSmall.RData")
model4 <- zipModelAlt


method <- c("iNim-RW", "MCMC")
ret <- rbind(c(model1$inlamcmc$mcmc.out$summary[,1], method[1]),
             #c(model2$inlamcmc$mcmc.out$summary[,1], method[2]),
             c(model3$mcmc$mcmc.out$summary[,1], method[2]))%>%
  as.data.frame()

retSD <- rbind(c(model1$inlamcmc$mcmc.out$summary[,3], method[1]),
               #c(model2$inlamcmc$mcmc.out$summary[,3], method[2]),
               c(model3$mcmc$mcmc.out$summary[,3], method[2]))%>%
  as.data.frame()

colnames(ret)[6] <- colnames(retSD)[6]  <- "method"
names <- colnames(ret)

#extract samples
method <- c("iNim2-RW")
ret1 <- c(model4[[1]]$mcmc.out$summary[,1], method)
retSD1 <- c(model4[[1]]$mcmc.out$summary[,3], method)
#colnames(ret1) <- colnames(retSD1) <- "method"
#save results
resMean <- cbind(rbind(ret,
                       ret1,
                       coefsEst),
                 "mean")
colnames(resMean)[7] <- "metrics"

resSD <- cbind(rbind(retSD,
                     retSD1,
                     coefsSD),
               "sd")

colnames(resSD) <- colnames(resMean)

allResults <- rbind(resMean,
                    resSD)%>%
  filter(!method %in% c("iNim-AFSS"))

# save results
write.csv(allResults, file = "results/zipResults.csv")

method <- c("iNim-RW", "MCMC")
ret <- rbind(cbind(model1$inlamcmc$mcmc.out$samples, method[1]),
             #cbind(model2$inlamcmc$mcmc.out$samples, method[2]),
             cbind(model3$mcmc$mcmc.out$samples, method[2]))%>%
  as.data.frame()
colnames(ret)[6] <- "method"

method <- c("iNim2-RW")
ret1 <- cbind(model4[[1]]$mcmc.out$samples, method)
colnames(ret1) <- colnames(ret)

columns_to_plot <- c("beta[1]",  "beta[2]", "beta[3]",
                     "gamma[1]", "gamma[2]")

allResults <- rbind(ret,
                    ret1)%>%
  as.data.frame()%>%
  mutate(method = as.factor(method))%>%
  mutate_at(vars(all_of(columns_to_plot)), as.numeric)

plotList <- list()

coefsEst[1:5] <- as.numeric(coefsEst[1:5])
for(i in seq_along(columns_to_plot)){

  col_greek <- parse(text = paste(columns_to_plot[i]))

  plot <- ggplot(allResults, aes(x = .data[[columns_to_plot[i]]], linetype = as.factor(method), col = as.factor(method)))+
    geom_density()+
    labs(title = col_greek)+
    xlab("")+
    theme_bw()+
    geom_vline(xintercept = as.numeric(coefsEst[i]))+
    scale_color_manual(values = fill.colors)+
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.5, "cm"))
  #print(plot)

  plotList[[i]] <- plot
}

zipPlot <- ggpubr::ggarrange(plotlist = plotList,
                                 ncol = 3,
                                 nrow = 2,
                                 common.legend = TRUE,
                                 legend = "bottom")


ggsave("results/binNmix.png",
       plot = binNmixPlot,
       width = 7,
       height = 6,
       units = "in")

##########
# Spatial Occupancy model
##########
load("spatialOccupancyModel/inlaWithNimble/occSpatialModelSmall1.RData")
model1 <- occSpatialModel

load("spatialOccupancyModel/inlaWithNimble/occSpatialModelSmall3.RData")
model3 <- occSpatialModel


load("spatialOccupancyModel/inlaWithNimbleAlt/spatialOccupancyAlt.RData")
model4 <- spatialModelAlt


load("spatialOccupancyModel/dataSimulated.RData")
coefsEst <- c(-1, 1, -0.45, 2, -2, 2.197, dataSimulated$variance.RF, dataSimulated$true_psi_fs, dataSimulated$theta.RF, "truth")

method <- c("iNim-RW", "MCMC")
ret <- rbind(c(model1$inlamcmc$mcmc.out$summary[1:7,2], mean(model1$inlamcmc$mcmc.out$summary[-c(1:7),2]), NA,method[1]),
             #c(model2$inlamcmc$mcmc.out$summary[,3], method[2]),
             c(model3$mcmc$mcmc.out$summary[1:7,2], mean(model3$mcmc$mcmc.out$summary[-c(1:7),2]),NA,method[2]))%>%
  as.data.frame()

retSD <- rbind(c(model1$inlamcmc$mcmc.out$summary[1:7,3], mean(model1$inlamcmc$mcmc.out$summary[-c(1:7),3]), NA, method[1]),
               #c(model2$inlamcmc$mcmc.out$summary[,3], method[2]),
               c(model3$mcmc$mcmc.out$summary[1:7,3], mean(model3$mcmc$mcmc.out$summary[-c(1:7),3]),NA, method[2]))%>%
  as.data.frame()

colnames(ret)[8:10] <- colnames(retSD)[8:10]  <- c("psi.fs","theta2","method")
names <- colnames(ret)


method <- c("iNim2-RW")
ret1 <- c(model4[[1]]$mcmc.out$summary[2:3,2],model4[[1]]$mcmc.out$summary[1,2], model4[[1]]$mcmc.out$summary[5:6,2],model4[[1]]$mcmc.out$summary[7,2],mean(model4[[1]]$mcmc.out$summary[-c(1:8),2]), model4[[1]]$mcmc.out$summary[8,2],method)
ret1SD <- c(model4[[1]]$mcmc.out$summary[2:3,3],model4[[1]]$mcmc.out$summary[1,3], model4[[1]]$mcmc.out$summary[5:6,3],model4[[1]]$mcmc.out$summary[7,3],mean(model4[[1]]$mcmc.out$summary[-c(1:8),3]), model4[[1]]$mcmc.out$summary[8,3],method)
names(ret1) <- names(ret1SD) <- names

resMean <- cbind(rbind(ret,
                       ret1,
                       coefsEst),
                 "mean")
colnames(resMean)[11] <- "metrics"

resSD <- cbind(rbind(retSD,
                     ret1SD),
               "sd")

colnames(resSD) <- colnames(resMean)

allResults <- rbind(resMean,
                    resSD)

# save results
write.csv(allResults, file = "results/spatialOccupancy.csv")
