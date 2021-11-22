runs <- read.csv("param_values.csv")

for (x in c("run", "slope", "slope_patches", "init_indivs", "init_patches", "num_gens")) {
    runs[[x]] <- NA
}
new_cols <- c("num_individuals", "num_juveniles", "num_patches", "num_hops_mean", "num_hops_sd", 
              "prop_dispersers", "distance_mean", "distance_sd", "children_mean", 
              "children_sd")
for (x in new_cols) runs[[x]] <- NA

for (j in 1:nrow(runs)) {
    logs <- list.files(as.character(runs$id[j]), "sim_.*.log", full.names=TRUE)
    if (length(logs) == 1) {
        runs$run[j] <- logs[1]
        x <- read.csv(logs[1], comment="#")
        for (cn in new_cols) {
            runs[[cn]][j] <- x[[cn]][nrow(x)]
        }
        runs$slope[j] <- coef(lm(num_individuals ~ generation, data=x))['generation']
        runs$slope_patches[j] <- coef(lm(num_patches ~ generation, data=x))['generation']
        runs$init_indivs[j] <- x$num_individuals[1]
        runs$init_patches[j] <- x$num_patches[1]
        runs$num_gens[j] <- x$generation[nrow(x)]
    } else if (length(logs) > 1) {
        stop(paste(c("More than one log file for", runs$id[j])))
    }
}

write.csv(runs, file="results.csv", row.names=FALSE)
