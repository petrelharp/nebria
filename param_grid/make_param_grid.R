library(jsonlite)
library(tidyverse)

default_params <- list(
    DEBUG = FALSE,
    NUM_GENS = 500,
    STEPSIZE = 1,
    POP_SIZE = 100,
    P_D = 0.2,
    YEAR_SHAPE = 1.5,
    DISPERSAL_SIGMA = 1.0
)

if (TRUE) {
    # Do first round of sims to get posterior
    default_params$NUM_GENS <- 500
    default_params$START_TIME_AGO <- default_params$NUM_GENS
    values <- list(
        POP_SIZE = c(40, 200), 
        DISPERSAL_SIGMA = c(0.2, 2),
        P_D = 2^-c(1, 6),
        YEAR_SHAPE = c(1, 2)
    )

    basedir <- "./post_500"
    nreps <- 50
    base_param_values <- do.call(expand.grid, values)
    param_values <- data.frame(lapply(base_param_values, function (x) {
                c(x, min(x) + runif(nreps - length(x)) * (max(x) - min(x)))
        }))
    param_values$id <- sprintf("run%06d", 1:nrow(param_values))
    dir.create(basedir, showWarnings=FALSE)
    write.csv(param_values, file=file.path(basedir, "param_values.csv"), row.names=FALSE)

}

if (FALSE) {
    # A few parameter values that give number of simulated patches close to 250
    default_params$NUM_GENS <- 21000
    default_params$START_TIME_AGO <- default_params$NUM_GENS
    basedir <- "./three_sims"
    param_values <- data.frame(POP_SIZE = c(69.7, 310.8, 111.9),
                               DISPERSAL_SIGMA = c(1.87, 0.28, 0.75),
                               P_D = c(0.48, 0.26, 0.08),
                               YEAR_SHAPE = c(1.95, 1.59, 1.26)
    )
}

if (FALSE) {
  # A few parameter values that give number of simulated patches close to 250
  default_params$NUM_GENS <- 21000
  default_params$START_TIME_AGO <- default_params$NUM_GENS
  basedir <- "./two_sims"
  param_values <- data.frame(POP_SIZE = c(57.10027, 87.82954),
                             DISPERSAL_SIGMA = c(0.5333791, 1.2358037),
                             P_D = c(0.44232404, 0.04112197),
                             YEAR_SHAPE = c(1.451843, 1.044548)
  )
}

if (FALSE) {
  set.seed(1003)
  # Draw 30 paramater values from the posterior distribution from the 500 simulations
  default_params$NUM_GENS <- 21000
  default_params$START_TIME_AGO <- default_params$NUM_GENS
  basedir <- "./post_21000"
  
  # Posterior samples
  post_500_res <- read.csv("post_500/posterior_samples.csv")
  
  param_values <- slice_sample(post_500_res, n = 30) %>% select(!X)
  param_values$id <- sprintf("run%06d", 1:nrow(param_values))

  dir.create(basedir, showWarnings=FALSE)
  write.csv(param_values, file=file.path(basedir, "param_values.csv"), row.names=FALSE)
}

if (FALSE) {
  seed <- 1004
  set.seed(seed)
  # Draw paramater values from the posterior distribution from the 500 simulations
  default_params$NUM_GENS <- 21000
  default_params$START_TIME_AGO <- default_params$NUM_GENS
  basedir <- "./post_21000b"
  
  # Posterior samples
  post_500_res <- read.csv("post_500/posterior_samples.csv") %>% rename(id=X)
  
  param_values <- slice_sample(post_500_res, n = 10)
  param_values$id <- sprintf("run%06d", param_values$id)

  dir.create(basedir, showWarnings=FALSE)
  write.csv(param_values, file=file.path(basedir, sprintf("param_values_%d.csv", seed)), row.names=FALSE)
}


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
        file.symlink(file.path("..", "..", "..", f), file.path(this_dir, f))
    }
}

