library(tidyverse)
source("../data/helpers.R")

# N/S split is on this longitude
long_break <- -118.68

# the data
observed_stats <- read_csv(file.path("../observed", "observe.stats.csv"))
observed_stats$site_name <- transl(observed_stats$site_name)
observed_pairstats <- read_csv(file.path("../observed", "observe.pairstats.csv"))
observed_pairstats$loc1 <- transl(observed_pairstats$loc1)
observed_pairstats$loc2 <- transl(observed_pairstats$loc2)
observed_pairstats$pair_name <- make_names(short_name(observed_pairstats$loc1), short_name(observed_pairstats$loc2))
sample_locs <- read_csv("../sample_locs.csv")
sample_locs$site_name <- transl(sample_locs$site_name)
stopifnot(nrow(observed_stats) == nrow(sample_locs))
observed_stats <- merge(observed_stats, sample_locs)
stopifnot(nrow(observed_stats) == nrow(sample_locs))
observed_stats$short_name <- short_name(observed_stats$site_name)
observed_stats$group <- ifelse(observed_stats$longitude < long_break, "N", "S")

observed_pairstats <- left_join(observed_pairstats, sample_locs, join_by(loc1 == site_name)) %>%
  left_join(sample_locs, join_by(loc2 == site_name), suffix = c("1", "2")) %>%
  select(loc1, loc2, dxy, pair_name, pct1, pct2) %>%
  mutate(pct_dist = abs(pct1 - pct2))

observed_pairstats$group <- paste(
  observed_stats$group[match(short_name(observed_pairstats$loc1), observed_stats$short_name)],
  observed_stats$group[match(short_name(observed_pairstats$loc2), observed_stats$short_name)],
  sep="-"
)


# Combine two simulation folders
#basedir1 = "post_21000"
#basedir2 = "post_21000_2022-12-05"

#all_stats1 <- read.csv(file.path(basedir1, "stats_all.csv"))
#all_stats2 <- read.csv(file.path(basedir2, "stats_all.csv"))

#all_pairstats1 <- read.csv(file.path(basedir1, "pairstats_all.csv"))
#all_pairstats2 <- read.csv(file.path(basedir2, "pairstats_all.csv"))

# Combine
#all_stats <- bind_rows(all_stats1, all_stats2)
#all_pairstats <- bind_rows(all_pairstats1, all_pairstats2)

basedir = "post_21000_2022-12-05"
all_stats <- read_csv(file.path(basedir, "stats_all.csv"))
all_pairstats <- read_csv(file.path(basedir, "pairstats_all.csv"))

# Split string specifying replicate into strings for sim, recap, mut, and sample
# Pattern to extract 

split_rep <- function(rep){
  match_pattern = "sim_([:digit:]+)_.+/.+_([:digit:]+)_([:digit:]+)_([:digit:]+)"
  sim  <- as.character(sapply(str_match_all(rep, match_pattern), "[", 2))
  recap  <- as.character(sapply(str_match_all(rep, match_pattern), "[", 3))
  mut  <- as.character(sapply(str_match_all(rep, match_pattern), "[", 4))
  sample  <- as.character(sapply(str_match_all(rep, match_pattern), "[", 5))
  split <- tibble(sim, recap, mut, sample)
  return(split)
}

all_stats <- all_stats %>% mutate(split_rep(rep),
                                  short_name = short_name(site_name))

rep_info <- data.frame( rep=unique(all_stats$rep) )
rep_info <- rep_info %>%  mutate(split_rep(rep))

for (pn in c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns", "mut_rate", "POP_SIZE", "P_D", "DISPERSAL_SIGMA",  "YEAR_SHAPE")) {
  rep_info[[pn]] <- NA
  for (j in 1:nrow(rep_info)) {
    this_rep <- rep_info$rep[j]
    x <- all_stats[[pn]][all_stats$rep == this_rep]
    stopifnot(length(unique(x)) == 1)
    rep_info[[pn]][j] <- x[1]
  }
}

# do the pairwise stats
all_pairstats <- all_pairstats %>% mutate(split_rep(rep),
                                        pair_name = make_names(short_name(loc1), short_name(loc2)))
# remove self-dxy
all_pairstats <- subset(all_pairstats, loc1 != loc2)

# make matrix of all reps x all stats
stats_names <- c(
  sprintf("het_%s", unique(short_name(all_stats$site_name))),
  sprintf("dxy_%s", unique(all_pairstats$pair_name))
)
stats <- matrix(NA, nrow=length(unique(all_stats$rep)), ncol=length(stats_names))
rownames(stats) <- rep_info$rep
colnames(stats) <- stats_names
for (j in 1:nrow(stats)) {
  this_rep <- rownames(stats)[j]
  x <- subset(all_stats, rep == this_rep)
  het_cols <- sprintf("het_%s", x$short_name)
  stats[j, het_cols] <- x$het
  y <- subset(all_pairstats, rep == this_rep)
  dxy_cols <- sprintf("dxy_%s", y$pair_name)
  stats[j, dxy_cols] <- y$dxy
}

observed <- rep(NA, ncol(stats))
names(observed) <- colnames(stats)
observed[sprintf("het_%s", observed_stats$short_name)] <- observed_stats$het
ut <- (observed_pairstats$loc1 != observed_pairstats$loc2)
observed[sprintf("dxy_%s", observed_pairstats$pair_name[ut])] <- observed_pairstats$dxy[ut]

# Write to file
results_folder <- file.path(basedir, "cleaned_results")
#dir.create(results_folder)
# Observed heterozygosity and pairwise stats
write.table(observed,  file.path(results_folder, "observed.csv"), sep = ",")
# Simulated heterozygosity and pairwise stats
write.csv(stats, file.path(results_folder, "stats.csv"))
# Parameters for each simulation and replicate
write.csv(rep_info, file.path(results_folder, "rep_info.csv"))
# Sampling locations
write.csv(sample_locs, file.path(results_folder, "sample_locs.csv"))
# Stats for all simulations
write.csv(all_stats, file.path(results_folder, "all_stats.csv"))
# Pairwise stats for all simulations
write.csv(all_pairstats, file.path(results_folder, "all_pairstats.csv"))

write.csv(observed_stats, file.path(results_folder, "observed_stats.csv"))
write.csv(observed_pairstats, file.path(results_folder, "observed_pairstats.csv"))
