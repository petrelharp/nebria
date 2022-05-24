library(randomForestSRC)
library(tidyverse)
set.seed(244)

all_stats <- read.csv("post_21000/stats_all.csv")
all_pairstats <- read.csv("post_21000/pairstats_all.csv")

all_stats$recap <- sapply(strsplit(all_stats$rep, "_"), "[", 5)
all_stats$mut <- sapply(strsplit(all_stats$rep, "_"), "[", 6)
all_stats$sample <- sapply(strsplit(all_stats$rep, "_"), "[", 7)
all_stats$sim <- sapply(strsplit(all_stats$rep, "_"), "[", 3)

all_pairstats$recap <- sapply(strsplit(all_pairstats$rep, "_"), "[", 5)
all_pairstats$mut <- sapply(strsplit(all_pairstats$rep, "_"), "[", 6)
all_pairstats$sample <- sapply(strsplit(all_pairstats$rep, "_"), "[", 7)
all_pairstats$sim <- sapply(strsplit(all_pairstats$rep, "_"), "[", 3)

one_sim <- filter(all_stats, sim == 2633463977134)
one_sim2 <- filter(all_pairstats, sim == 2633463977134)


# 27 sites
# Choose(27, 2) + 27 pairwise comparisons of sites
# 24 SLiM simulations
# Per simulation 5 recap reps
# Per recap rep 5 mutation reps
# Per mutation rep 5 sample reps
# 24*5*5*5 = 3000

cols_to_merge <- c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns", "POP_SIZE", "P_D", "YEAR_SHAPE", "sim", "sample", "recap", "mut")
obs_stats <- "het"
obs_pairstats <- "dxy"

# Make one row for each of the 3000 sim/recap/mut/sample replicates
wide_stats <- pivot_wider(all_stats, id_cols = all_of(cols_to_merge),
            names_from = site_name,
            names_glue = "{site_name}_{.value}",
            values_from = all_of(obs_stats))
wide_pairstats <- pivot_wider(all_pairstats, id_cols = all_of(cols_to_merge),
                              names_from = c("loc1", "loc2"),
                              names_glue = "{loc1}_{loc2}_{.value}",
                              values_from = all_of(obs_pairstats))
wide_all <- inner_join(wide_stats, wide_pairstats, by = all_of(cols_to_merge))

# Random forest
sim_rep <- sample(wide_all$sim, 1)
train <- filter(wide_all, sim == sim_rep)
test <- filter(wide_all, sim != sim_rep)

#params_to_predict <- c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns", "POP_SIZE", "P_D", "YEAR_SHAPE")
params_to_predict <- c("POP_SIZE", "P_D", "YEAR_SHAPE")
stats_allsites <- c(grep("het", names(wide_all), value = TRUE), grep("dxy", names(wide_all), value = TRUE))

rf_data <- select(train, all_of(c(params_to_predict, stats_allsites)))
model <- rfsrc(Multivar(POP_SIZE, P_D, YEAR_SHAPE) ~ . , rf_data, forest = TRUE)

prediction <- predict(model$forest, test)
predicted_params <- prediction$regrOutput
plot(test$POP_SIZE, predicted_params$POP_SIZE$predicted)
