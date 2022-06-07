library(tidyverse)

basedir = "post_21000/"
all_stats <- read.csv(file.path(basedir, "stats_all.csv"))

all_stats$sim <- factor(sapply(strsplit(all_stats$rep, "_"), "[", 3))
all_stats$recap <- factor(sapply(strsplit(all_stats$rep, "_"), "[", 5))
all_stats$mut <- factor(sapply(strsplit(all_stats$rep, "_"), "[", 6))
all_stats$sample <- sapply(strsplit(all_stats$rep, "_"), "[", 7)

one_sim <- filter(all_stats, sim == all_stats$sim[76572])
ggplot(one_sim, aes(x = mut, y = het, color = sample)) +
  geom_point() +
  facet_wrap(~recap)

# Look at variation within a few simulations
pdf(file="within-sim-variation.pdf", width=8, height=10, pointsize=10)
layout(matrix(1:6, ncol=2))
site_list <- c("Conness Lake", "Sixty Lakes", "Army Pass")
for (sim_id in c("2633463977134", "1843169324096", "1745513977159", "3424803977142", "2633493977134", "2607102226138")) {
    plot(0, type='n', xlim=c(0, 125), ylim=range(0, all_stats$het), xlab='Index', ylab='heterozygosity',
         main=paste("sim ", sim_id)
    )
    for (k in seq_along(site_list)) {
        site <- site_list[k]
        with(droplevels(subset(all_stats, sim == sim_id & site_name == site)), {
             points(het,
                    pch=20,
                    col=k,
            )
            abline(v=0.5 + which(diff(as.numeric(recap)) > 0), lty=1)
            abline(v=0.5 + which(diff(as.numeric(mut)) > 0), lty=3)
        })
    }
    legend("bottomright",
           pch=c(NA, NA, rep(20, length(site_list))),
           lty=c(1, 3, rep(NA, length(site_list))),
           col=c(1, 1, seq_along(site_list)),
           legend=c("recap reps", "mutation reps", site_list),
           bg='white'
    )
}
dev.off()

# do the pairwise stats
all_pairstats <- read.csv(file.path(basedir, "pairstats_all.csv"))
all_pairstats$sim <- factor(sapply(strsplit(all_pairstats$rep, "_"), "[", 3))
all_pairstats$recap <- factor(sapply(strsplit(all_pairstats$rep, "_"), "[", 5))
all_pairstats$mut <- factor(sapply(strsplit(all_pairstats$rep, "_"), "[", 6))
all_pairstats$sample <- sapply(strsplit(all_pairstats$rep, "_"), "[", 7)

