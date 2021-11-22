runs <- read.csv("results.csv")
pairs(runs[,c("POP_SIZE", "DISPERSAL_SIGMA", "P_D", "YEAR_SHAPE", "num_individuals", "slope", "num_gens")])

varnames <- c("POP_SIZE", "DISPERSAL_SIGMA", "P_D", "YEAR_SHAPE")
cor(runs[,varnames], runs$slope)

sigma_vals <- cut(runs$DISPERSAL_SIGMA, breaks=7)

plot(slope / num_individuals ~ jitter(POP_SIZE), data=runs, col=as.numeric(sigma_vals), pch=20)
legend('topleft', col=1:nlevels(sigma_vals), pch=20, pt.cex=2,
       legend=sprintf("dispersal sigma = %0.2f", levels(sigma_vals)))


plot(num_individuals ~ jitter(POP_SIZE), data=runs, col=match(DISPERSAL_SIGMA, sigma_vals), pch=20)
legend('topleft', col=seq_along(sigma_vals), pch=20, pt.cex=2, legend=sprintf("dispersal sigma = %0.2f", sigma_vals))

