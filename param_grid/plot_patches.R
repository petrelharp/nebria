library(tidyverse)

# Want a plot of num_individuals and num_patches over generations for every simulation

raster_switch_times_ago <- c(21, 
                        (17 + 14.7)/2,
                        (14.7 + 12.9)/2,
                        (12.9 + 11.7)/2,
                        (11.7 + 8.326)/2,
                        (8.326 + 4.2)/2,
                        (4.2 + 0.3)/2,
                        0) * 1000
raster_switch_gens <- 21000 - raster_switch_times_ago
basedir = "two_sims/"

sim_results <- read.csv(paste0(basedir, "results_all_gens.csv"))

ggplot(sim_results, aes(x = generation, y = num_individuals)) +
  geom_point() +
  facet_wrap(~filename) +
  geom_vline(xintercept = raster_switch_gens)
ggplot(sim_results, aes(x = generation, y = num_patches)) +
  geom_line() +
  facet_wrap(~filename) +
  geom_vline(xintercept = raster_switch_gens)

