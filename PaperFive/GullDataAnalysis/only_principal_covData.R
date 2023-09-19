source("nimble.R")

only_principal_cov <- run_simulations(dataout, "fixed-covariate", TRUE)

save(only_principal_cov , file="estimated_data_only_principal_cov2.RData")
