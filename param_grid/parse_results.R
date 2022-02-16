basedir = "two_sims/"

runs <- read.csv(paste0(basedir, "param_values.csv"))

for (x in c("run", "slope", "slope_patches", "init_indivs", "init_patches", "num_gens")) {
    runs[[x]] <- NA
}

# new_cols <- c("num_individuals", "num_juveniles", "num_patches", "num_hops_mean", "num_hops_sd",
#                "prop_dispersers", "distance_mean", "distance_sd", "children_mean",
#                "children_sd")
new_cols <- c("num_individuals", "num_juveniles", "num_patches")
for (x in new_cols) runs[[x]] <- NA

for (j in 1:nrow(runs)) {
    logs <- list.files(paste0(basedir, runs$id[j]), "sim_.*.log", full.names=TRUE)
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

write.csv(runs, file=paste0(basedir, "results.csv"), row.names=FALSE)

# Combine all simulation log files into one
read_plus <- function(flnm, ...){
  read_csv(flnm, show_col_types = FALSE, ...) %>% mutate(filename = flnm)
}

run_folders <- paste0(basedir, runs$id)
log_files <- list.files(run_folders, "sim_.*.log", recursive = TRUE, full.names=TRUE)
all_gens <- log_files %>% map_df(~read_plus(., comment = "#"))

write.csv(all_gens, file=paste0(basedir, "results_all_gens.csv"), row.names=FALSE)

