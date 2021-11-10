library(jsonlite)

default_params  <- fromJSON("params.json")
default_params$MAX_SIZE <- 1e6
default_params$NUM_GENS <- 10

values <- list(
    POP_SIZE = 10 * 2^(2:5),
    DISPERSAL_SIGMA = seq(0.2, 2, length.out=5),
    P_D = 2^-(2:6),
    YEAR_SHAPE = c(1, 1.5, 2)
)

param_values <- do.call(expand.grid, values)
param_values$id <- sprintf("run%06d", 1:nrow(param_values))

basedir <- "param_grid"
dir.create(basedir, showWarnings=FALSE)
write.csv(param_values, file=file.path(basedir, "param_values.csv"), row.names=FALSE)
writeLines(toJSON(default_params, pretty=TRUE), file.path(basedir, "default_params.json"))

setup_files <- c("BA13_15.png", "current.png", "EH8_12.png", "geo_only_suitability.png", "HS15_17.png", "LH0_4.png", "MH4_8.png", "LGM17_21.png")

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
