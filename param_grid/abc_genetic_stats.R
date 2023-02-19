library(tidyverse)
library(abc)

# Observed heterozygosity and pairwise stats
observed <- read.csv("post_21000_2022-12-05/cleaned_results/observed.csv")
# Simulated heterozygosity and pairwise stats
stats <- read.csv("post_21000_2022-12-05/cleaned_results/stats.csv")
# Parameters for each simulation and replicate
rep_info <- read.csv("post_21000_2022-12-05/cleaned_results/rep_info.csv")
rep_info <- rep_info |> mutate(across(c(sim, recap, mut, sample), factor))
# Sampling locations
sample_locs <- read.csv("post_21000_2022-12-05/cleaned_results/sample_locs.csv")

# Check that observed and simulated summary statistics are in the same order
all(rownames(observed) == names(dplyr::select(stats, !X)))
# Check that info and stats are in the same order
all(rep_info$rep == stats$X)

# Parameters we want to estimate
est_params <- c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns", "mut_rate", "POP_SIZE", "P_D", "YEAR_SHAPE")
abc_obs <- observed/mean(observed$x)
abc_sims <- dplyr::select(stats, !X)/rowMeans(dplyr::select(stats, !X))
abc_res2 <- abc(abc_obs, dplyr::select(rep_info, est_params), abc_sims , tol = 0.1, method = "loclinear") 
abc_res <- abc(observed, dplyr::select(rep_info, est_params), dplyr::select(stats, !X) , tol = 0.1, method = "loclinear") 

#!!!! Why are there NA mutation rates? !!!!
est_params <- c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns", "mut_rate", "POP_SIZE", "P_D", "YEAR_SHAPE")
na_reps <- which(is.na(rep_info$mut_rate))
abc_obs <- observed$x
abc_sims <- dplyr::select(stats, !X) %>% filter(!is.na(rep_info$mut_rate))
abc_rep_info <- dplyr::select(rep_info, est_params) %>% filter(!is.na(mut_rate))
abc_res <- abc(abc_obs, abc_rep_info, abc_sims , tol = 0.1, method = "rejection") 
#abc_res <- abc(abc_obs, abc_rep_info, abc_sims , tol = 0.1, method = "neuralnet") 

# Rejection method results

rej_post <- data.frame(abc_res$unadj.values) %>% pivot_longer(everything(), names_to = "parameter")
priors <- rep_info %>% dplyr::select(est_params) %>% pivot_longer(everything(), names_to = "parameter")

#adj_post <- data.frame(abc_res$adj.values) %>% pivot_longer(everything(), names_to = "parameter")

ggplot(rej_post, aes(x = value, y = after_stat(density))) +
  geom_histogram() +
  geom_histogram(data = priors, fill = "blue", alpha = 0.5) +
  facet_wrap(~parameter, scales = "free")

ggplot(rej_post, aes(x = value, y = after_stat(density))) +
  geom_density() +
  geom_density(data = priors, fill = "blue", alpha = 0.5) +
  facet_wrap(~parameter, scales = "free")

pairs(abc_res$unadj.values)

rej_res <- data.frame(abc_res$unadj.values)
ggplot(rej_res, aes(x = P_D)) +
  geom_histogram(bins = 100) +
  geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)
ggplot(rej_res, aes(x = POP_SIZE, bins = 100)) +
  geom_histogram() + 
  geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)
ggplot(rej_res, aes(x = YEAR_SHAPE)) +
  geom_histogram() +
  geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)
ggplot(rej_res, aes(x = NE)) +
  geom_histogram() +
  geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)
ggplot(rej_res, aes(x = NE)) +
  geom_histogram() +
  geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)

for(est_p in est_params){
  ggplot(rej_res, aes(x = est_p)) +
    geom_histogram() +
    geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)
}



# ggplot(adj_post, aes(x = value)) +
#   geom_histogram() +
#   facet_wrap(~parameter)

adj_res <- data.frame(abc_res$adj.values)
rej_res <- data.frame(abc_res$unadj.values)

# Save results for adjusted posterior

write.csv(rej_res, file = paste0(basedir, "posterior_samples.csv"))


ggplot(adj_res, aes(x = P_D)) +
  geom_histogram(bins = 100) +
  geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)
ggplot(adj_res, aes(x = POP_SIZE, bins = 100)) +
  geom_histogram() + 
  geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)
ggplot(adj_res, aes(x = YEAR_SHAPE)) +
  geom_histogram() +
  geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)
ggplot(adj_res, aes(x = NE)) +
  geom_histogram() +
  geom_histogram(data = rep_info, fill = "blue", alpha = 0.5)

pairs(adj_res)


# Hold out a simulation and do ABC
# Simulation we're using as observation
set.seed(100)
(sim_observed_num <- sample(levels(rep_info$sim), 1))
(sim_observed_rows <- which(rep_info == sim_observed_num))
(sim_observed <- stats[sim_observed_rows,])
