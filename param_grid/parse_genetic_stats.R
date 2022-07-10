library(jsonlite)
library(tidyverse)

basedir = "post_21000/"
nsites = 27

stats_files <- list.files(basedir, pattern = "*\\.stats.csv", recursive = TRUE, full.names = TRUE)
base_files <- gsub("\\.stats.csv", "", stats_files)

merge_stats <- function(base_file, ending){
  # ending is .stats.csv or .pairstats.csv
    json_file <- paste0(base_file, ".json")
    stopifnot(file.exists(json_file))
  recap_params <- fromJSON(json_file)
  run_dir <- grep("run*", strsplit(base_file, "/")[[1]], value = TRUE)
  sim_params <- fromJSON(paste0(basedir, run_dir,"/","params.json"))
  stats <- read.csv(paste0(base_file, ending))
  merged <- merge(stats, data.frame(recap_params)) %>% 
    merge(., data.frame(sim_params)) %>%
    mutate(rep = base_file)
  return(merged)
}
all_stats <- base_files %>% map_df(~merge_stats(., ending = ".stats.csv"))

if(nrow(all_stats) != length(base_files)*nsites){
  print("Error: wrong number of rows")
}

#all_stats$recap_rep <- sapply(strsplit(all_stats$rep, "_"), "[", 4)
#all_stats$mut_rep <- sapply(strsplit(all_stats$rep, "_"), "[", 5)
#all_stats$sample_rep <- sapply(strsplit(all_stats$rep, "_"), "[", 6)
#all_stats$sim_rep <- sapply(strsplit(all_stats$rep, "_"), "[", 3)

write.csv(all_stats, file=paste0(basedir, "stats_all.csv"), row.names=FALSE)

pairstats_files <- list.files(basedir, pattern = "*\\.pairstats.csv", recursive = TRUE, full.names = TRUE)
base_files <- gsub("\\.pairstats.csv", "", pairstats_files)

all_pairstats <- base_files %>% map_df(~merge_stats(., ending = ".pairstats.csv"))
if(nrow(all_pairstats) != length(base_files)*(choose(nsites, 2) + nsites)){
  print("Error: wrong number of rows")
}

write.csv(all_pairstats, file=paste0(basedir, "pairstats_all.csv"), row.names=FALSE)
