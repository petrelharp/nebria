library(tidyverse)
library(abc)
library(jsonlite)

source("../data/helpers.R")

basedir <- "post_21000_2022-12-05"
# Observed heterozygosity and pairwise stats
observed_df <- read.csv(file.path(basedir, "cleaned_results/observed.csv"))
observed <- observed_df$x
names(observed) <- rownames(observed_df)
# Simulated heterozygosity and pairwise stats
stats <- read.csv(file.path(basedir,"cleaned_results/stats.csv"))
stats <- dplyr::select(stats, !X)
# Parameters for each simulation and replicate
rep_info <- read.csv(file.path(basedir,"cleaned_results/rep_info.csv"))
rep_info <- rep_info |> mutate(across(c(sim, recap, mut, sample), factor))

# Check that observed and simulated summary statistics are in the same order
all(rownames(observed) == names(stats))
# Check that info and stats are in the same order
all(rep_info$rep == stats$X)

#!!!! Why are there NA mutation rates? !!!!
est_params <- c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns", "mut_rate", "POP_SIZE", "P_D", "DISPERSAL_SIGMA", "YEAR_SHAPE")
na_reps <- which(is.na(rep_info$mut_rate))
abc_obs <- observed
abc_sims <- stats %>% filter(!is.na(rep_info$mut_rate))
abc_rep_info <- dplyr::select(rep_info, est_params) %>% filter(!is.na(mut_rate))
#abc_res2 <- abc(abc_obs, abc_rep_info, abc_sims , tol = 0.1, method = "rejection")
# This only works with a higher tolerance, why?
abc_res <- abc(abc_obs, abc_rep_info, abc_sims , tol = 0.2, method = "loclinear")
#abc_res3 <- abc(abc_obs, abc_rep_info, abc_sims , tol = 0.1, method = "ridge")

rej_post <- data.frame(abc_res$unadj.values) %>% pivot_longer(everything(), names_to = "parameter")
adj_post <- data.frame(abc_res$adj.values) %>% pivot_longer(everything(), names_to = "parameter")
priors <- rep_info %>% dplyr::select(est_params) %>% pivot_longer(everything(), names_to = "parameter")

estimates <- apply(abc_res$adj.values, 2, median)

# Setup a new simulation with parameters drawn from the posterior distribution

set.seed(1000)
nsamples <- 10
draws <- sample_n(data.frame(abc_res$adj.values), nsamples)

# Remove draws with negative parameters
valid_draws <- filter(draws, if_all(names(draws), ~ . > 0))
print(paste0(nrow(draws) - nrow(valid_draws), " row(s) removed because of negative param values"))

default_params <- list(
  DEBUG = FALSE,
  NUM_GENS = 500,
  STEPSIZE = 1,
  POP_SIZE = 100,
  P_D = 0.2,
  YEAR_SHAPE = 1.5,
  DISPERSAL_SIGMA = 1.0
)

default_params$NUM_GENS <- 21000
default_params$START_TIME_AGO <- default_params$NUM_GENS
datestring = format(Sys.time(), "%Y-%m-%d")
basedir <- paste0("./estimated_sims_", datestring)

param_values <- valid_draws
param_values$id <- sprintf("run_%s_%06d", datestring, (1:nrow(param_values)))

dir.create(basedir, showWarnings=FALSE)
setup_files <- c("geo_layers")
write.csv(valid_draws, file.path(basedir,"posterior_draws.csv"))

for(j in 1:nrow(param_values)) {
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

#######Sean code
#ABC:Possible method values are "rejection", "loclinear", "neuralnet" and "ridge".

#First conduct Cross-validation for  parameter inference
cv.resR <- cv4abc(param=abc_rep_info, abc_sims, nval=5,tols=c(0.05,0.10,0.15,0.2,0.25,0.3,0.35), method="rejection")
cv.resL <- cv4abc(param=abc_rep_info, abc_sims, nval=5,tols=c(.15,0.2,0.25), method="loclinear")
cv.resRd <- cv4abc(param=abc_rep_info, abc_sims, nval=5,tols=c(.15,0.2,0.25), method="ridge")
df <- as.data.frame(summary(cv.resR))
df <- rbind(df,summary(cv.resL))
df <- rbind(df,summary(cv.resRd))

#which.min(df[,1]) #use this to examine tolerance level with lowest error for each variable
#For rejection and loclinear, tol 0.25. Ridge does not perform well.
#Do a neural net (but cv4abc function doesn't work well) at tol=0.25

abc_resR <- abc(abc_obs, abc_rep_info, abc_sims , tol = 0.25, method = "rejection")
abc_resL <- abc(abc_obs, abc_rep_info, abc_sims , tol = 0.25, method = "loclinear")
abc_resN <- abc(abc_obs, abc_rep_info, abc_sims , tol = 0.25, method = "neuralnet",sizenet=2)
#had to reduce the sizenet to 2 (from defaulty 5) due to weights error

summary(abc_resR)
summary(abc_resL)
summary(abc_resN)

#Save parameters from the rejection and neural net fit.
rej_tbl <- summary(abc_resR)
write.table(rej_tbl,file="abc_results/parameters_rejection.txt",sep="\t")

nnet_tbl <- summary(abc_resN)
write.table(nnet_tbl,file="abc_results/parameters_neuralnet.txt",sep="\t")

posterior <- data.frame(abc_resN$adj.values) %>% rename("Average population size per patch" = POP_SIZE,
                                                        "Probability of dispersal" = P_D,
                                                        "Dispersal distance" = DISPERSAL_SIGMA,
                                                        "Yearly variation in habitat quality" = YEAR_SHAPE,
                                                        "Mutation rate" = mut_rate,
                                                        T2=T2, T1=T1, CS=CS, AS=AS, NE=NE, Na=Na, Nc=Nc, Ns=Ns) %>%
  pivot_longer(everything(), names_to = "Parameter")

estimates_vec <- apply(abc_resN$adj.values, 2, median)
names(estimates_vec) <- c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns", 
                         "Mutation rate", "Average population size per patch", 
                         "Probability of dispersal",
                         "Dispersal distance", 
                         "Yearly variation in habitat quality")
estimates <- data.frame(Parameter = names(estimates_vec), value = estimates_vec)

prior <- abc_rep_info %>% rename("Average population size per patch" = POP_SIZE,
                                 "Probability of dispersal" = P_D,
                                 "Dispersal distance" = DISPERSAL_SIGMA,
                                 "Yearly variation in habitat quality" = YEAR_SHAPE,
                                 "Mutation rate" = mut_rate,
                                 T2=T2, T1=T1, CS=CS, AS=AS, NE=NE, Na=Na, Nc=Nc, Ns=Ns) %>%
  pivot_longer(everything(), names_to = "Parameter" )
msprime_params <- c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns")
slim_params <- c("Average population size per patch", "Probability of dispersal",
                 "Dispersal distance", "Yearly variation in habitat quality")
mut_params <- c("Mutation rate")
posterior_msprime <- filter(posterior, Parameter %in% msprime_params)
prior_msprime <- filter(prior, Parameter %in% msprime_params)
est_msprime <- filter(estimates, Parameter %in% msprime_params)
posterior_slim <- filter(posterior, Parameter %in% slim_params)
prior_slim <- filter(prior, Parameter %in% slim_params)
est_slim <- filter(estimates, Parameter %in% slim_params)
posterior_mut <- filter(posterior, Parameter %in% mut_params)
prior_mut <- filter(prior, Parameter %in% mut_params)
est_mut <- filter(estimates, Parameter %in% mut_params)

ggplot(prior, aes(x = value)) +
  geom_density() +
  geom_density(data = posterior, aes(x = value), color = "blue") +
  geom_vline(aes(xintercept = value), data = estimates, color = "blue") +
  facet_wrap(~Parameter, scales="free")

dens_cols <- c("prior" = "black", "posterior" = "blue")
ggplot(prior_msprime, aes(x = value, fill = "prior")) +
  geom_density(alpha = 0.5) +
  geom_density(data = posterior_msprime, aes(x = value, fill = "posterior"), alpha = 0.5) +
  geom_vline(aes(xintercept = value), data = est_msprime, color = "blue") +
  facet_wrap(~Parameter, scales="free", ncol = 2) +
  scale_fill_manual(name = "", values = dens_cols) +
  theme(text=element_text(size=18),
        axis.text.x = element_text(angle = 90))
ggsave("../generate_figures/msprime_params.png")
ggplot(prior_slim, aes(x = value, fill = "prior")) +
  geom_density(alpha = 0.5) +
  geom_density(data = posterior_slim, aes(x = value, fill = "posterior"), alpha = 0.5) +
  geom_vline(aes(xintercept = value), data = est_slim, color = "blue") +
  facet_wrap(~Parameter, scales="free") +
  scale_fill_manual(name = "", values = dens_cols) +
  theme(text=element_text(size=18),
        axis.text.x = element_text(angle = 90))
ggsave("../generate_figures/slim_params.png")
ggplot(prior_mut, aes(x = value, fill = "prior")) +
  geom_density(alpha = 0.5) +
  geom_density(data = posterior_mut, aes(x = value, fill = "posterior"), alpha = 0.5) +
  geom_vline(aes(xintercept = value), data = est_mut, color = "blue") +
  facet_wrap(~Parameter, scales="free") +
  scale_fill_manual(name = "", values = dens_cols) +
  theme(text=element_text(size=18),
        axis.text.x = element_text(angle = 90))
ggsave("../generate_figures/mut_params.png")

