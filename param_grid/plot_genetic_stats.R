library(tidyverse)

basedir = "post_21000/"
all_stats <- read.csv(paste0(basedir, "stats_all.csv"))
all_pairstats <- read.csv(paste0(basedir, "pairstats_all.csv"))

ggplot(all_stats, aes(x = rep, color = rep, y = het)) +
  geom_point()
