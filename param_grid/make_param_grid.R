library(jsonlite)

default_params  <- fromJSON("params.json")
default_params$MAX_SIZE <- 1e6
default_params$NUM_GENS <- 5
values <- list(
    POP_SIZE = c(40, 320), 
    DISPERSAL_SIGMA = c(0.2, 2),
    P_D = 2^-c(1, 6),
    YEAR_SHAPE = c(1, 2)
)

nreps <- 300
base_param_values <- do.call(expand.grid, values)
param_values <- data.frame(lapply(base_param_values, function (x) {
            c(x, min(x) + runif(nreps - length(x)) * (max(x) - min(x)))
    }))
param_values$id <- sprintf("run%06d", 1:nrow(param_values))

basedir <- "."
dir.create(basedir, showWarnings=FALSE)
write.csv(param_values, file=file.path(basedir, "param_values.csv"), row.names=FALSE)
writeLines(toJSON(default_params, pretty=TRUE), file.path(basedir, "default_params.json"))

setup_files <- c("geo_layers")

for (j in 1:nrow(param_values)) {
    this_dir <- file.path(basedir, param_values$id[j])
    dir.create(this_dir, recursive=TRUE, showWarnings=FALSE)
    params <- default_params
    for (xn in names(param_values)) {
        if (xn %in% names(params)) {
            params[[xn]] <- param_values[[xn]][j]
        }
    }
    writeLines(toJSON(params, pretty=TRUE), file.path(this_dir, "params.json"))
    for (f in setup_files) {
        file.symlink(file.path("..", "..", f), file.path(this_dir, f))
    }
}

