runs <- read.csv("param_values.csv")

for (x in c("run", "slop", "init_indivs")) {
    runs[[x]] <- NA
}
new_cols <- c("num_individuals", "num_juveniles", "num_hops_mean", "num_hops_sd", 
              "prop_dispersers", "distance_mean", "distance_sd", "children_mean", 
              "children_sd")
for (x in new_cols) runs[[x]] <- NA

for (j in 1:nrow(runs)) {
    logs <- list.files(as.character(runs$id[j]), "sim_.*.log", full.names=TRUE)
    runs$run[j] <- logs[1]
    x <- read.csv(logs[1], comment="#")
    for (cn in new_cols) {
        runs[[cn]][j] <- x[[cn]][nrow(x)]
    }
    runs$slope[j] <- coef(lm(num_individuals ~ generation, data=x))['generation']
    runs$init_indivs[j] <- x$num_individuals[1]
}

write.csv(runs, file="results.csv", row.names=FALSE)
