# This script is used to put together the results from the CV for the insect
# groups.
# It requires the loading of data from the crossValidation folder and simple 
# manipulation of the results

# library needed to run this script
library(dplyr)

# Load data
load("crossValidation/bumblebees/estimateCrossValidate.RData")
bb <- estimatesCaseStudy
load("crossValidation/hoverflies/estimateCrossValidate.RData")
hv <- estimatesCaseStudy
load("crossValidation/solitarybees/estimateCrossValidate.RData")
sb <- estimatesCaseStudy

#put all data together
allResults <- c(bb, 
                hv, 
                sb)


#Extract CV values

cvVals <- lapply(allResults, function(x){
  ret <- x$CVvalue
  diff <- ret[1] - ret[2]
  retVal <- c(ret[1:2], diff)
})%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(InsectGroup = c(rep("bumblebees", 5),
                         rep("hoverflies", 5),
                         rep("solitarybees", 5)),
         Model = rep(c("GCMSH", "GCMCO", "SOM", "IDMSH", "IDMCO"), 3))

colnames(cvVals)[1:3] <- c('CVvalue', 'GC_IDM', 'SO_IDM')

write.csv(cvVals, file = "Figures/cvVals.csv", row.names = FALSE)



# Marginal contribution
# I do this by hand

#Bumblebees
fitCount <- cvVals[4:5,3] - cvVals[3,1]
PanTrap <- cvVals[4:5,2] - cvVals[1:2,1]

# Hoverflies
fitCount <- cvVals[9:10,3] - cvVals[8,1]
PanTrap <- cvVals[9:10,2] - cvVals[6:7,1]


# Solitarybees
fitCount <- cvVals[14:15,3] - cvVals[13,1]
PanTrap <- cvVals[14:15,2] - cvVals[11:12,1]
