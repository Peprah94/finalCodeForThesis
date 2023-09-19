source("nimble.R")

main <- run_simulations(dataout, "main", TRUE)

save(main, file="estimated_data_main2.RData")