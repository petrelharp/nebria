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
observed_pairstats$pct_dist <- abs(sample_locs$pct[observed_pairstats$loc1]
                                   - sample_locs$pct[observed_pairstats$loc2])
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
all_stats$sim <- factor(sapply(str_match_all(all_stats$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 2))
all_stats$recap <- factor(sapply(str_match_all(all_stats$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 3))
all_stats$mut <- factor(sapply(str_match_all(all_stats$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 4))
all_stats$sample <- factor(sapply(str_match_all(all_stats$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 5))
all_stats$site_name <- transl(all_stats$site_name)
all_stats$short_name <- short_name(all_stats$site_name)

rep_info <- data.frame( rep=unique(all_stats$rep) )
rep_info$sim <- factor(sapply(str_match_all(rep_info$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 2))
rep_info$recap <- factor(sapply(str_match_all(rep_info$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 3))
rep_info$mut <- factor(sapply(str_match_all(rep_info$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 4))
rep_info$sample <- factor(sapply(str_match_all(rep_info$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 5))

for (pn in c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns", "mut_rate", "POP_SIZE", "P_D", "YEAR_SHAPE")) {
  rep_info[[pn]] <- NA
  for (j in 1:nrow(rep_info)) {
    this_rep <- rep_info$rep[j]
    x <- all_stats[[pn]][all_stats$rep == this_rep]
    stopifnot(length(unique(x)) == 1)
    rep_info[[pn]][j] <- x[1]
  }
}

# do the pairwise stats
all_pairstats$sim <- factor(sapply(str_match_all(all_pairstats$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 2))
all_pairstats$recap <- factor(sapply(str_match_all(all_pairstats$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 3))
all_pairstats$mut <- factor(sapply(str_match_all(all_pairstats$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 4))
all_pairstats$sample <- factor(sapply(str_match_all(all_pairstats$rep, "sim_([:digit:]+)_stats/stats_([:digit:]+)_([:digit:]+)_([:digit:]+)"), "[", 5))
all_pairstats$pair_name <- make_names(short_name(all_pairstats$loc1), short_name(all_pairstats$loc2))

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
dir.create(results_folder)
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



