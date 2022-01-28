runs <- read.csv("sim_runs/results.csv")
varnames <- c("POP_SIZE", "DISPERSAL_SIGMA", "P_D", "YEAR_SHAPE")

pairs(runs[,c(varnames, "num_individuals", "num_patches", "slope", "slope_patches", "num_gens")],
      pch=20, col=ifelse(is.na(runs$slope), 'red', 
                         ifelse(runs$slope < 8000, 'black', 'green'))
)

cor(runs[,varnames], runs$slope_patches, use='pairwise')

sigma_vals <- cut(runs$DISPERSAL_SIGMA, breaks=7)

plot(slope / num_individuals ~ jitter(POP_SIZE), data=runs, col=as.numeric(sigma_vals), pch=20)
legend('topleft', col=1:nlevels(sigma_vals), pch=20, pt.cex=2,
       legend=sprintf("dispersal sigma = %s", levels(sigma_vals)))


plot(num_individuals ~ jitter(POP_SIZE), data=runs, col=sigma_vals, pch=20)
legend('topleft', col=1:nlevels(sigma_vals), pch=20, pt.cex=2,
       legend=sprintf("dispersal sigma = %s", levels(sigma_vals)))


the_lm <- lm(slope ~ POP_SIZE + DISPERSAL_SIGMA + P_D + YEAR_SHAPE, data=runs)

# number of juveniles per adult
with(runs, hist(num_juveniles/num_individuals))

# number of individuals per patch
with(runs, hist(num_individuals/num_patches, breaks=30))

