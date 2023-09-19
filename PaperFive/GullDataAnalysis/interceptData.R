source("nimble.R")

intercept <- run_simulations(dataout, "intercept", TRUE)

save(intercept, file="estimated_data_intercept2.RData")