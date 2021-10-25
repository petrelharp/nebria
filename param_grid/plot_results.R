runs <- read.csv("results.csv")
pairs(runs[,c("POP_SIZE", "DISPERSAL_SIGMA", "P_D", "YEAR_SHAPE", "num_individuals", "slope")])

