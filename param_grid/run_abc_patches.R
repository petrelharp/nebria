library(abc)
library(tidyverse)

runs <- read.csv("sim_runs_new_maps/results.csv")

# Want a posterior distribution for, all with uniform priors over the range
# P_D (1/2 to 1/64)
# DISPERSAL_SIGMA (0.2 to 2.0)
# POP_SIZE, (40 to 320)
# YEAR_SHAPE (1 to 2)

# Simulated parameters
sim_params <- select(runs, P_D, DISPERSAL_SIGMA, POP_SIZE, YEAR_SHAPE)
# Simulated summary statistics (number of occupied patches)
# Add normal noise with standard deviation 100
sim_patches <- pmax(0, runs$num_patches + rnorm(nrow(runs), 0, sd = 100))

# Observed summary statistic (number of occupied patches)
obs_patches <- 250
abc_res <- abc(target = obs_patches, param = sim_params, sumstat = sim_patches,
               tol = 0.1, method = "loclinear")

# Rejection method results

rej_post <- data.frame(abc_res$unadj.values) %>% pivot_longer(everything(), names_to = "parameter")
adj_post <- data.frame(abc_res$adj.values) %>% pivot_longer(everything(), names_to = "parameter")

ggplot(rej_post, aes(x = value)) +
  geom_histogram() +
  facet_wrap(~parameter)
ggplot(adj_post, aes(x = value)) +
  geom_histogram() +
  facet_wrap(~parameter)

adj_res <- data.frame(abc_res$adj.values)

# Save results for adjusted posterior

write.csv(adj_res, file = "sim_runs_new_maps/posterior_samples.csv")

ggplot(adj_res, aes(x = DISPERSAL_SIGMA)) +
  geom_histogram()
ggplot(adj_res, aes(x = P_D)) +
  geom_histogram()
ggplot(adj_res, aes(x = POP_SIZE)) +
  geom_histogram()
ggplot(adj_res, aes(x = YEAR_SHAPE)) +
  geom_histogram()

varnames <- c("POP_SIZE", "DISPERSAL_SIGMA", "P_D", "YEAR_SHAPE")
pairs(runs[,varnames],
      pch=ifelse(runs$num_gens == 40, 20, 1), col=ifelse(runs$num_patches < 500 & runs$num_patches > 100, 'red', 'black'),
      )
plot(runs$POP_SIZE, runs$DISPERSAL_SIGMA, col=ifelse(runs$num_patches < 500 & runs$num_patches > 100, 'red', 'black'))
#identify(runs$POP_SIZE, runs$DISPERSAL_SIGMA)

ggplot(runs, aes(x = POP_SIZE, y = DISPERSAL_SIGMA, color = abs(num_patches - 250))) +
  geom_point()

close_patches <- filter(runs, num_patches < 450 & num_patches > 150)
ggplot(close_patches, aes(x = num_individuals)) +
  geom_histogram()
ggplot(close_patches, aes(x = num_juveniles)) +
  geom_histogram()

# Parameters that give observed patches around 250

params_to_check <- filter(runs, num_patches <= 260 & num_patches >= 240) %>% select(POP_SIZE, DISPERSAL_SIGMA, P_D, YEAR_SHAPE, id, num_patches)
