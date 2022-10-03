library(tidyverse)
source("../data/helpers.R")

# N/S split is on this longitude
long_break <- -118.68

# the data
observed_stats <- read.csv(file.path("../observed", "observe.stats.csv"))
observed_stats$site_name <- transl(observed_stats$site_name)
observed_pairstats <- read.csv(file.path("../observed", "observe.pairstats.csv"))
observed_pairstats$loc1 <- transl(observed_pairstats$loc1)
observed_pairstats$loc2 <- transl(observed_pairstats$loc2)
observed_pairstats$pair_name <- make_names(short_name(observed_pairstats$loc1), short_name(observed_pairstats$loc2))
sample_locs <- read.csv("../sample_locs.csv")
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

basedir = "post_21000/"
all_stats <- read.csv(file.path(basedir, "stats_all.csv"))

all_stats$sim <- factor(sapply(strsplit(all_stats$rep, "_"), "[", 3))
all_stats$recap <- factor(sapply(strsplit(all_stats$rep, "_"), "[", 5))
all_stats$mut <- factor(sapply(strsplit(all_stats$rep, "_"), "[", 6))
all_stats$sample <- sapply(strsplit(all_stats$rep, "_"), "[", 7)
all_stats$site_name <- transl(all_stats$site_name)
all_stats$short_name <- short_name(all_stats$site_name)

rep_info <- data.frame( rep=unique(all_stats$rep) )
rep_info$sim <- factor(sapply(strsplit(rep_info$rep, "_"), "[", 3))
rep_info$recap <- factor(sapply(strsplit(rep_info$rep, "_"), "[", 5))
rep_info$mut <- factor(sapply(strsplit(rep_info$rep, "_"), "[", 6))
rep_info$sample <- sapply(strsplit(rep_info$rep, "_"), "[", 7)

for (pn in c("T2", "T1", "CS", "AS", "NE", "Na", "Nc", "Ns", "POP_SIZE", "P_D", "YEAR_SHAPE")) {
    rep_info[[pn]] <- NA
    for (j in 1:nrow(rep_info)) {
        this_rep <- rep_info$rep[j]
        x <- all_stats[[pn]][all_stats$rep == this_rep]
        stopifnot(length(unique(x)) == 1)
        rep_info[[pn]][j] <- x[1]
    }
}

if (FALSE) {
    # look at one sim
    one_sim <- filter(all_stats, sim == all_stats$sim[76572])
    ggplot(one_sim, aes(x = mut, y = het, color = sample)) +
      geom_point() +
      facet_wrap(~recap)
}

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

# make a good order for plotting things:
# first hets in order of pct,
# then dxys in order of pct distance
# within/between groups split in the middle
key <- rep(NA, length(observed))
names(key) <- names(observed)
key[sprintf("het_%s", observed_stats$short_name)] <- observed_stats$pct - 1e10
ut <- (observed_pairstats$loc1 != observed_pairstats$loc2)
key[sprintf("dxy_%s", observed_pairstats$pair_name[ut])] <- (
    1e8 * observed_stats$pct[short_name(observed_pairstats$loc1[ut])]
    + observed_stats$pct[short_name(observed_pairstats$loc2[ut])]
)
plot_ord <- rank(key)
stopifnot(length(plot_ord) == length(observed))

# normalization
norm <- function (x) {
    ut <- grepl("dxy", names(x))
    stopifnot(sum(ut) == choose(27, 2))
    xm <- mean(x[ut])
    return(x/xm)
}

stopifnot(all(names(observed) == colnames(stats)))
## boxplots
plot_ranges <- function (stats, observed, plot_ord, norm=identity, ...) {
    stats <- do.call(rbind, lapply(1:nrow(stats), function (j) norm(stats[j,])))
    observed <- norm(observed)
    statsum <- apply(stats, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
    hets <- grepl("het", names(observed))
    # plot_ord <- rank(observed - 1000 * grepl("het", colnames(stats)))
    par(mar=c(10, 3, 3, 1))
    plot(0, type='n', xlim=range(plot_ord), ylim=range(statsum, observed), xlab='', xaxt='n', ylab='normalized heterozygosity', ...)
    segments(x0=plot_ord, y0=statsum["2.5%",], y1=statsum["97.5%",])
    segments(x0=plot_ord, y0=statsum["25%",], y1=statsum["75%",], lwd=3)
    points(x=plot_ord, y=statsum["50%",], col='black', cex=3, pch=20)
    points(x=plot_ord, y=observed, col='blue', cex=3, pch=20)
    axis(1, at=plot_ord, labels=names(observed), las=3, cex.axis=0.7)
    loc1_labels <- tapply(plot_ord,
                          gsub("_.*$", "", gsub("dxy_", "", names(observed))),
                          mean
    )
    loc1_labels <- loc1_labels[!grepl("het", names(loc1_labels))]
    text(loc1_labels, 0.4, labels=names(loc1_labels), pos=1)
}

pdf(file="normalized_stats.pdf", width=36, height=6, pointsize=8)
plot_ranges(stats, observed, plot_ord, main=basedir, norm=norm)
dev.off()

pdf(file="unnormalized_stats.pdf", width=36, height=6, pointsize=8)
plot_ranges(stats, observed, plot_ord, main=basedir)
dev.off()

if (FALSE) {
    # TODO: adjust this code to make plots for each location n-s
    plot_ranges2 <- function (stats, pair_stats, ...) {
        # from ../plot_stats.R
        stats <- stats[order(stats$pct),]
        yscale <- 1.0 * max(pair_stats$dxy, na.rm=TRUE)
        par(mar=c(5, 6, 3, 5)+.1)
        plot(0, type='n',
             xlim=range(stats$pct),
             xlab='PCT distance [m]',
             ylim=c(0.75, nrow(stats)-0.75) * yscale, yaxt='n',
             ylab='',
             main="dxy"
        )
        axis(2, las=2, tick=FALSE,
             at=(1:nrow(stats) - 0.5) * yscale,
             labels=stats$short_name
        )
        abline(h=(0:nrow(stats)) * yscale, lwd=2, col=adjustcolor('black', 0.5))
        axis(4, las=2, tick=TRUE,
             at=(0:nrow(stats)) * yscale,
             labels=sprintf(
                "%0.1e",
                (0:nrow(stats)) * yscale - 2 * rep(0:nrow(stats), each=2)[1:(nrow(stats)+1)] * yscale
            )
        )
        for (k in 1:nrow(stats)) {
            ref = stats$short_name[k]
            ii = match(make_names(ref, stats$short_name), pair_stats$pair_name)
            y = yscale * (k-1) + pair_stats$dxy[ii]
            abline(h=y[k], lty=3, col=adjustcolor(k, 0.5))
            lines(stats$pct, y, col=k, lty=k, type='b', pch=20)
            points(stats$pct[k], y[k], pch="*", cex=5, col=k)
        }
    }

    pdf(file="normalized_stats2.pdf", width=6.5, height=15, pointsize=10)
        plot_ranges2(observed_stats, observed_pairstats)
    dev.off()
}

## PCA

# Remove row 247 because it has NA for the pairwise stats. I don't know why, I need to figure it out
# Actually replace row 247 with a duplicate of row 246 because the rest of the code depends on having the same number of rows

stats[247,] <- stats[246,]
sum(is.na(stats)) # should be 0
pc_stats <- prcomp(stats)
for (k in 1:4) {
    rep_info[[paste0("PC", k)]] <- pc_stats$x[,k]
}

# Convert observed data into PC space
pc_observed <- predict(pc_stats, data.frame(t(observed)))

pdf("sims_pca.pdf", width=12, height=5, pointsize=10)
    ggplot(rep_info, aes(x=PC1, y=PC2, col=P_D)) + 
        geom_point() + facet_wrap(~cut(POP_SIZE,3))
dev.off()

ggplot(rep_info) +
  geom_point(aes(x = PC1, y = PC2, col = sim)) +
  geom_point(aes(x = pc_observed[,'PC1'], y = pc_observed[,'PC2']))
ggplot(rep_info) +
  geom_point(aes(x = PC3, y = PC4, col = sim)) +
  geom_point(aes(x = pc_observed[,'PC3'], y = pc_observed[,'PC4']))

ggplot(rep_info) +
  geom_point(aes(x = PC1, y = PC2, col = recap)) +
  geom_point(aes(x = pc_observed[,'PC1'], y = pc_observed[,'PC2'])) +
  facet_wrap(~sim)
ggplot(rep_info) +
  geom_point(aes(x = PC3, y = PC4, col = recap)) +
  geom_point(aes(x = pc_observed[,'PC3'], y = pc_observed[,'PC4'])) +
  facet_wrap(~sim)

plot(pc_stats$sdev^2/sum(pc_stats$sdev^2))
one_sim <- filter(rep_info, sim == 3696507884660)
nrow(one_sim)
length(unique(one_sim$mut))
length(unique(one_sim$rep))
length(unique(one_sim$recap))

