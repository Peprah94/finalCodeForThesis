source("nimble.R")

variable <- run_simulations(dataout, "variable", TRUE)

save(variable, file="estimated_data_variable2.RData")
