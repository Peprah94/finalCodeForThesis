#This script further format the gull dataset with the ML prediction scores

library(readr)
library(dplyr)

#Load the prediction scores
observations <- read_csv("GullDataAnalysis/observations.csv")%>%
  filter(`Prediction 1` %in% c("Larus argentatus",
                               "Larus marinus",
                               "Larus canus", 
                               "Larus fuscus"))

rr <- cbind(observations$`Prediction 1 score`,
            observations$`Prediction 2 score`,
            observations$`Prediction 3 score`,
            observations$`Prediction 4 score`,
            observations$`Prediction 5 score`)%>%
  as.matrix()


data_obs <- data.frame(gbifID = rep(observations$gbifID, 5),
                       scientificName = rep(observations$scientificName, 5),
                       occurenceID = rep(observations$occurrenceID, 5),
                       predicted_species = c(observations$`Prediction 1`,
                                             observations$`Prediction 2`,
                                             observations$`Prediction 3`,
                                             observations$`Prediction 4`,
                                             observations$`Prediction 5`
                       ),
                       predicted_score = c(observations$`Prediction 1 score`,
                                           observations$`Prediction 2 score`,
                                           observations$`Prediction 3 score`,
                                           observations$`Prediction 4 score`,
                                           observations$`Prediction 5 score`))%>%
  tidyr::pivot_wider(., id_cols = c("gbifID", "occurenceID"),
                     names_from = predicted_species,
                     values_from = predicted_score,
                     values_fill = 0
  ) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(others = sum(c_across(-c(1:6)), na.rm = T))%>%
  dplyr::ungroup()%>%
  #dplyr::select(1:5,8, 377)%>%
  dplyr::select(1:6, 320)%>%
  dplyr::arrange(gbifID)

colnames(data_obs) <- c("gbifID", "occurenceID", "herring",
                        "great black_backed",
                        "common",  
                        "lesser black_backed", 
                        "other")

#reorder column names
col_order <- c("gbifID", "occurenceID", "common",
               "great black_backed",
               "herring",
                "lesser black_backed", 
               "other")
data_obs <- data_obs[, col_order]


#Load gull dataset
load("gull_df.RData")

#Add the ML data
datas <- gull_df %>%
  dplyr::inner_join(data_obs, by = c("gbifID"))

# Format covariates
gull <- datas%>%
  mutate(experience = as.factor(ifelse(previous_obs + 1 > 9, 10, previous_obs + 1)),
         altitude = (as.numeric(scale(alt))),
         temperature = (as.numeric(scale(bio1))),
         precipitation = (as.numeric(scale(bio12)))#,
         #scale_dist = as.numeric(log10(distance^2))
  )%>%
  filter_at(vars(altitude, temperature, precipitation), all_vars(!is.na(.)))%>%
  filter_at(vars(altitude, temperature, precipitation), all_vars(!is.infinite(.)))

#PLotting covariates for exploration
gull%>%
  dplyr::select(experience, altitude, temperature, precipitation, 
                truth, reported)%>%
  GGally::ggpairs()

# Extracting validation and training dataset
validation_data_with_distance <- gull%>%
  filter(year == 2022 & countryCode != "NO")

training_data_with_distance <- gull%>%
  filter(year!= 2022 | countryCode == "NO")

trainingDataTable <- table(training_data_with_distance$truth,
      training_data_with_distance$reported)

write.csv(trainingDataTable,"trainingDataSummaryTable.csv", row.names = FALSE)

validationDataTable <- table(validation_data_with_distance$truth,
                           validation_data_with_distance$reported)
write.csv(validationDataTable,"validationDataSummaryTable.csv", row.names = FALSE)

Y <- as.numeric(as.factor(c(training_data_with_distance$reported, 
                            validation_data_with_distance$reported)))

C <- as.numeric(as.factor(c(training_data_with_distance$truth, 
                            validation_data_with_distance$truth)))

index_for_splitting <- seq(nrow(training_data_with_distance)+1, 
                           nrow(training_data_with_distance)+nrow(validation_data_with_distance))

Validation_C <- as.numeric(as.factor(validation_data_with_distance$truth))
Validation_Y <- as.numeric(as.factor(validation_data_with_distance$reported))

validation_indices_for_mismatch = index_for_splitting[which(Validation_C != Validation_Y)]
validation_indices_for_match = index_for_splitting[which(Validation_C == Validation_Y)]
validation_mismatchC = C[validation_indices_for_mismatch]
validation_mismatchY = Y[ validation_indices_for_mismatch]
validation_matchC = C[ validation_indices_for_match]
validation_matchY = Y[validation_indices_for_match]

#table(Validation_C,Validation_Y)
C[index_for_splitting] <- NA
#Y[1:n.visit, index_for_splitting] <- NA

omega <- rbind(training_data_with_distance[,300:304], 
               validation_data_with_distance[,300:304])%>%
  as.matrix()

#covariates
altitude <- c(training_data_with_distance$altitude, validation_data_with_distance$altitude)
distance <- c(training_data_with_distance$scale_dist, validation_data_with_distance$scale_dist)
precipitation <- c(training_data_with_distance$precipitation, validation_data_with_distance$precipitation)
mean_temperature <- c(training_data_with_distance$temperature, validation_data_with_distance$temperature)

cov = rbind(altitude, 
            precipitation)





#unique user
users <- as.numeric(as.factor( c(training_data_with_distance$recordedBy, validation_data_with_distance$recordedBy)))

#experience
#experience <- c(training_data_with_distance$previous_obs, validation_data_with_distance$previous_obs)
experience <- c(training_data_with_distance$experience, validation_data_with_distance$experience)
no_previous_obs <- c(training_data_with_distance$previous_obs, validation_data_with_distance$previous_obs)


C <- matrix(as.numeric(as.factor(C)), nrow=1) # Unverified citizen science observation
Y <- matrix(as.numeric(as.factor(Y)), nrow=1) #Verified Citizen science observation



dataout <- list(C=C, #Verified data
                Y=Y, # Reported data
                omega=omega, # Confusion matrix
                cov=cov, # covariate
                cov_omega = experience,
                Validation_C = Validation_C,
                Validation_Y = Validation_Y,
                validation_indices_for_mismatch = validation_indices_for_mismatch,
                index_for_splitting = index_for_splitting,
                validation_mismatchC = validation_mismatchC,
                validation_mismatchY = validation_mismatchY,
                validation_indices_for_match = validation_indices_for_match,
                validation_matchC = validation_matchC,
                validation_matchY = validation_matchY
)

distance_to_road <- list(training_data_with_distance, 
                         validation_data_with_distance)

#load("distance_to_road.RData")
save(distance_to_road, file="distance_to_road_ML.RData")
save(dataout, file = "nimble_data_ML.RData")



  