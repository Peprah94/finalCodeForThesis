source("nimble.R")

ml_estimates <- run_simulations(dataout, "ML", TRUE)

save(ml_estimates , file="estimated_data_ml_estimates2.RData")
