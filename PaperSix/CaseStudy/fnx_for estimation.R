mbias <- function(true, estimate){
  ret <- mean(true-estimate) 
  return(ret)
}


mse <- function(true, estimate){
  ret <- mean((true-estimate)^2)
  return(ret)
}


hill_index <- function(q, pis){
    hill <- rowSums((pis^q))^(1/(1-q)) 
  return(hill)
}

#hill_index(2,pis)
#estimation of shannon index
shan_index <- function(pis){
  #shan <- exp(-rowSums(pis * log(pis)))
 shan <-  -rowSums(log(pis)*(pis), na.rm = TRUE)
  return(shan)
}

#Estimation of eveness
eveness <- function(pis){
  shan <- shan_index(pis)
  richness <- log(hill_index(0,pis))
  even <- shan/richness
return(even)
  }



mbias <- function(true, estimate){
  ret <- mean(true-estimate) 
  return(ret)
}


mse <- function(true, estimate){
  ret <- mean((true-estimate)^2)
  return(ret)
}

#Trials
#shan_index(pis)
#hill_index(2, pis)
#eveness(pis)