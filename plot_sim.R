#!/usr/bin/env R

args = commandArgs(TRUE)

usage = "
Usage:
    plot_sims.R (basename)

Will read in the following files:
 - sample_locs.csv : in this directory, at least needs
    * site_name
    * longitude
    * latitude
 - (basename).stats.csv : should be a CSV with the following columns (and possibly others):
    * site_name:
    * het: heterozygosity
 - (basename).pairstats.csv: a CSV with columns:
    * loc1: a name in 'site_name' from the stats file
    * loc2: a name in 'site_name' from the stats file
    * dxy: divergence

For background maps, you'll also need the `data/` subdirectory.
"

if (length(args) != 1) {
    stop(usage)
}

basename = gsub(".trees$", "", args[1])
outfile <- sprintf("%s.dxy_by_pct.png", basename)

library(png)

source("data/helpers.R")

sample_locs = read.csv("sample_locs.csv")
sample_locs$site_name <- transl(sample_locs$site_name)
sample_locs$short_name = gsub("[^a-z].*", "", tolower(sapply(strsplit(sample_locs$site_name, " "), "[", 1)))
rownames(sample_locs) = sample_locs$short_name
sample_locs <- sample_locs[order(sample_locs$pct),]

stats_data = read.csv(sprintf("%s.stats.csv", basename))
stats_data$site_name <- transl(stats_data$site_name)
stats_data = merge(sample_locs, stats_data, all=TRUE)

pairs_data = read.csv(sprintf("%s.pairstats.csv", basename))
pairs_data$loc1 <- transl(pairs_data$loc1)
pairs_data$loc2 <- transl(pairs_data$loc2)
pairs_data = pairs_data[,setdiff(colnames(pairs_data), "X")]

rownames(pairs_data) <- make_names(
                    sample_locs$short_name[match(pairs_data$loc1, sample_locs$site_name)],
                    sample_locs$short_name[match(pairs_data$loc2, sample_locs$site_name)]
               )

png(file=outfile, width=6.5, height=15, pointsize=10)
    yscale <- 1.0 * max(pairs_data$dxy, na.rm=TRUE)
    par(mar=c(5, 6, 3, 5)+.1)
    plot(0, type='n',
         xlim=range(sample_locs$pct),
         xlab='PCT distance [m]',
         ylim=c(0.75, nrow(sample_locs)-0.75) * yscale, yaxt='n',
         ylab='',
         main="dxy"
    )
    axis(2, las=2, tick=FALSE,
         at=(1:nrow(sample_locs) - 0.5) * yscale,
         labels=rownames(sample_locs)
    )
    abline(h=(0:nrow(sample_locs)) * yscale, lwd=2, col=adjustcolor('black', 0.5))
    axis(4, las=2, tick=TRUE,
         at=(0:nrow(sample_locs)) * yscale,
         labels=sprintf(
            "%0.1e",
            (0:nrow(sample_locs)) * yscale - 2 * rep(0:nrow(sample_locs), each=2)[1:(nrow(sample_locs)+1)] * yscale
        )
    )
    for (k in 1:nrow(sample_locs)) {
        ref = rownames(sample_locs)[k]
        ii = match(make_names(ref, rownames(sample_locs)), rownames(pairs_data))
        y = yscale * (k-1) + pairs_data$dxy[ii]
        abline(h=y[k], lty=3, col=adjustcolor(k, 0.5))
        lines(sample_locs$pct, y, col=k, lty=k, type='b', pch=20)
        points(sample_locs$pct[k], y[k], pch="*", cex=5, col=k)
    }
dev.off()

cat(sprintf("Wrote file to %s\n", outfile))

