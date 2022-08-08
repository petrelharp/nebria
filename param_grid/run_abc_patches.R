library(abc)
library(tidyverse)

basedir = "post_500/"

runs <- read.csv(paste0(basedir, "results.csv"))

# Want a posterior distribution for, all with uniform priors over the range
# P_D (1/2 to 1/64)
# DISPERSAL_SIGMA (0.2 to 2.0)
# POP_SIZE, (40 to 200)
# YEAR_SHAPE (1 to 2)

# Simulated parameters
sim_params <- dplyr::select(runs, P_D, DISPERSAL_SIGMA, POP_SIZE, YEAR_SHAPE)
# Simulated summary statistics (number of occupied patches)
# Add normal noise with standard deviation 100
sim_patches <- pmax(0, runs$num_patches + rnorm(nrow(runs), 0, sd = 0))

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
rej_res <- data.frame(abc_res$unadj.values)

# Save results for adjusted posterior

write.csv(rej_res, file = paste0(basedir, "posterior_samples.csv"))

ggplot(adj_res, aes(x = DISPERSAL_SIGMA)) +
  geom_histogram()
ggplot(adj_res, aes(x = P_D)) +
  geom_histogram()
ggplot(adj_res, aes(x = POP_SIZE)) +
  geom_histogram()
ggplot(adj_res, aes(x = YEAR_SHAPE)) +
  geom_histogram()

ggplot(adj_res, aes(x = DISPERSAL_SIGMA, y = POP_SIZE)) +
  geom_point()

pairs(adj_res)
pairs(rej_res)

