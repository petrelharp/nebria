#!/usr/bin/env Rscript
# Set up simulations based on posterior draws

args <- commandArgs(TRUE)
if (length(args) != 2) {
    stop("Usage: setup_sims.R <start row> <end row>")
}
start_row <- args[1]
end_row <- args[2]

library(jsonlite)

datestring = format(Sys.time(), "%Y-%m-%d")

param_values <- read.csv("abc_results/posterior_draws.csv", row.names=1)
param_values$id <- sprintf("run_%s_%06d", datestring, (1:nrow(param_values)))
param_values <- param_values[start_row:end_row,]

if (start_row > end_row || end_row > nrow(param_values)) {
    stop("Invalid start and/or end row")
}

default_params <- list(
  DEBUG = FALSE,
  NUM_GENS = 21000,
  START_TIME_AGO = 21000,
  STEPSIZE = 1
)

basedir <- paste0("./nn_sims_", datestring)


dir.create(basedir, showWarnings=FALSE)
setup_files <- c("geo_layers")
write.csv(param_values, file.path(basedir,"posterior_draws.csv"))

for(j in 1:nrow(param_values)) {
  this_dir <- file.path(basedir, param_values$id[j])
  dir.create(this_dir, recursive=TRUE, showWarnings=FALSE)
  ### WRONG!!!
  params <- default_params
  for (xn in names(param_values)) {
    params[[xn]] <- param_values[[xn]][j]
  }
  writeLines(jsonlite::toJSON(params, pretty=TRUE), file.path(this_dir, "params.json"))
  for (f in setup_files) {
    file.symlink(file.path("..", "..", "..", f), file.path(this_dir, f))
  }
}

