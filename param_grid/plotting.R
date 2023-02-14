library(tidyverse)
source("../data/helpers.R")
basedir <- "post_21000_2022-12-05"
# Observed heterozygosity and pairwise stats
observed_df <- read.csv(file.path(basedir, "cleaned_results/observed.csv"))
observed <- observed_df$x
names(observed) <- rownames(observed_df)
# Simulated heterozygosity and pairwise stats
stats <- read.csv(file.path(basedir,"cleaned_results/stats.csv")) %>% select(!X)
# Parameters for each simulation and replicate
rep_info <- read.csv(file.path(basedir,"cleaned_results/rep_info.csv"))
rep_info <- rep_info |> mutate(across(c(sim, recap, mut, sample), factor))
# Sampling locations
sample_locs <- read.csv(file.path(basedir,"cleaned_results/sample_locs.csv"))
# Stats for all simulations
all_stats <- read.csv(file.path(basedir,"cleaned_results/all_stats.csv"))
# Pairwise stats for all simulations
all_pairstats <- read.csv(file.path(basedir,"cleaned_results/all_pairstats.csv"))

observed_stats <- read.csv(file.path(basedir,"cleaned_results/observed_stats.csv"))
observed_pairstats <- read.csv(file.path(basedir,"cleaned_results/observed_pairstats.csv"))

all_stats$sim <- factor(all_stats$sim)
all_pairstats$sim <- factor(all_pairstats$sim)
subset_sims <- sample(all_stats$sim, 10)
small_all_stats <- filter(all_stats, sim %in% subset_sims)
# Plot values of observed statistics for each simulation replicate
ggplot() +
  geom_point(data = small_all_stats, aes(x = rep, y = het, col = sim)) +
  geom_point(data = observed_stats, aes(x = "observed", y = het))
filter(rep_info, sim %in% subset_sims) %>% select(mut_rate, sim)
small_pairstats <- filter(all_pairstats, sim %in% subset_sims)
ggplot() +
  geom_point(data = small_pairstats, aes(x = rep, y = dxy, col = sim)) +
  geom_point(data = observed_pairstats, aes(x = "observed", y = dxy))



## PCA

# Remove row 247 for post_21000 because it has NA for the pairwise stats. I don't know why, I need to figure it out
# Actually replace row 247 with a duplicate of row 246 because the rest of the code depends on having the same number of rows

#stats[247,] <- stats[246,]
sum(is.na(stats)) # should be 0

# Explore ways to normalize because first PC is negative and constant
pc_stats <- prcomp(stats[,-1], scale = TRUE)
for (k in 1:4) {
  rep_info[[paste0("PC", k)]] <- pc_stats$x[,k]
}

# Convert observed data into PC space
pc_observed <- predict(pc_stats, data.frame(t(observed)))

pdf("sims_pca.pdf", width=12, height=5, pointsize=10)
ggplot(rep_info, aes(x=PC1, y=PC2, col=P_D)) + 
  geom_point() + facet_wrap(~cut(POP_SIZE,3))
dev.off()

ggplot(rep_info) +
  geom_point(aes(x = PC1, y = PC2, col = sim)) +
  geom_point(aes(x = pc_observed[,'PC1'], y = pc_observed[,'PC2']))
ggplot(rep_info) +
  geom_point(aes(x = PC3, y = PC4, col = sim)) +
  geom_point(aes(x = pc_observed[,'PC3'], y = pc_observed[,'PC4']))

ggplot(rep_info) +
  geom_point(aes(x = PC1, y = PC2, col = recap)) +
  geom_point(aes(x = pc_observed[,'PC1'], y = pc_observed[,'PC2'])) +
  facet_wrap(~sim)
ggplot(rep_info) +
  geom_point(aes(x = PC3, y = PC4, col = recap)) +
  geom_point(aes(x = pc_observed[,'PC3'], y = pc_observed[,'PC4'])) +
  facet_wrap(~sim)

plot(pc_stats$sdev^2/sum(pc_stats$sdev^2))
one_sim <- filter(rep_info, sim == 3696507884660)
nrow(one_sim)
length(unique(one_sim$mut))
length(unique(one_sim$rep))
length(unique(one_sim$recap))

