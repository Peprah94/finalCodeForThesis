source("nimble.R")

op_cov_finter <- run_simulations(dataout, "fixec-intercov", TRUE)

save(op_cov_finter, file="estimated_data_op_cov_finter2.RData")
