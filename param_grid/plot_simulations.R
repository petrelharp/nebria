library(abc)
library(tidyverse)

basedir = "post_500/"

runs <- read.csv(paste0(basedir, "results.csv"))
runs <- filter(runs, !is.na(num_patches))

varnames <- c("POP_SIZE", "DISPERSAL_SIGMA", "P_D", "YEAR_SHAPE")
pairs(runs[,varnames],
      pch=ifelse(runs$num_gens == 40, 20, 1), col=ifelse(runs$num_patches < 500 & runs$num_patches > 100, 'red', 'black'),
)
plot(runs$POP_SIZE, runs$DISPERSAL_SIGMA, col=ifelse(runs$num_patches < 500 & runs$num_patches > 100, 'red', 'black'))
#identify(runs$POP_SIZE, runs$DISPERSAL_SIGMA)

ggplot(runs, aes(x = POP_SIZE, y = DISPERSAL_SIGMA, color = abs(num_patches - 250))) +
  geom_point()

close_patches <- filter(runs, num_patches < 450 & num_patches > 150)
ggplot(close_patches, aes(x = num_individuals)) +
  geom_histogram()
ggplot(close_patches, aes(x = num_juveniles)) +
  geom_histogram()

# Parameters that give observed patches around 250

params_to_check <- filter(runs, num_patches <= 400 & num_patches >= 100) %>% select(POP_SIZE, DISPERSAL_SIGMA, P_D, YEAR_SHAPE, id, num_patches)

# Filter to exclude runs with very large population sizes ( > 100000)
ggplot(runs,aes(x = num_individuals)) +
  geom_histogram()
sum(runs$num_individuals > 100000)
ggplot(runs, aes(x= num_individuals,y = num_patches )) +
  geom_point()
runs_smaller <- filter(runs, num_individuals < 100000)
ggplot(runs_smaller, aes(x = num_individuals, y = num_patches)) +
  geom_point()
ggplot(runs, aes(x = POP_SIZE, y = num_patches)) +
  geom_point()
ggplot(runs, aes(x = DISPERSAL_SIGMA, y = num_patches)) +
  geom_point()
ggplot(runs, aes(x = P_D, y = num_patches)) +
  geom_point()
ggplot(runs, aes(x = YEAR_SHAPE, y = num_patches)) +
  geom_point()

# Seems like population size should be between 0 and 200
ggplot(runs_smaller, aes(x = POP_SIZE, y = DISPERSAL_SIGMA, col=(num_patches < 500 & num_patches > 100))) +
  geom_point(show.legend = FALSE)

# Get a few parameter combinations with patches close to 250
filter(runs, num_patches < 270 & num_patches > 230)
