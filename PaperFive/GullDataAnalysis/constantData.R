source("nimble.R")

constant <- run_simulations(dataout, "constant", TRUE)

save(constant, file="estimated_data_constant2.RData")