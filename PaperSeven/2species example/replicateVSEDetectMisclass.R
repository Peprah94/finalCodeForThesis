
for(i in 1:10){
load(paste0("VSEDetectMisclassProbs/VSEDetectMisclass", i,".RData"))

save(mcmc.out, file = paste0("VSEDetectMisclass/VSEDetectMisclass",i+90,".RData"))
}
