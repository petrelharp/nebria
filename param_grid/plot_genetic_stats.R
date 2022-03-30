library(tidyverse)

basedir = "post_21000/"
all_stats <- read.csv(paste0(basedir, "stats_all.csv"))
#all_pairstats <- read.csv(paste0(basedir, "pairstats_all.csv"))

all_stats$recap_rep <- sapply(strsplit(all_stats$rep, "_"), "[", 4)
all_stats$mut_rep <- sapply(strsplit(all_stats$rep, "_"), "[", 5)
all_stats$sample_rep <- sapply(strsplit(all_stats$rep, "_"), "[", 6)
all_stats$sim_rep <- sapply(strsplit(all_stats$rep, "_"), "[", 3)

one_sim <- filter(all_stats, sim_rep == all_stats$sim_rep[76572])
ggplot(one_sim, aes(x = mut_rep, y = het, color = sample_rep)) +
  geom_point() +
  facet_wrap(~recap_rep)
