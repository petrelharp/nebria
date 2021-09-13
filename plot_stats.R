#!/usr/bin/env R

args = commandArgs(TRUE)

usage = "
Usage:
    plot_stats.R (basename)

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

library(png)
library(tidyverse)
library(sf)
library(stars)

source("data/mapping-fns.R", chdir=TRUE)

sample_locs = read.csv("sample_locs.csv")

stats_data = read.csv(sprintf("%s.stats.csv", basename))
stats_data = merge(sample_locs, stats_data, all=TRUE)

pairs_data = read.csv(sprintf("%s.pairstats.csv", basename))
pairs_data = pairs_data[,setdiff(colnames(pairs_data), "X")]

stats = st_as_sf(
            stats_data, agr="identity",
            coords = c("longitude", "latitude"), dim="XY",
            crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
)
stats$short_name = gsub("[^a-z].*", "", tolower(sapply(strsplit(stats$site_name, " "), "[", 1)))
rownames(stats) = stats$short_name

l_lines = lapply(1:nrow(pairs_data), function (k) {
                    i = match(pairs_data$loc1[k], stats$site_name)
                    j = match(pairs_data$loc2[k], stats$site_name)
                    stats[c(i,j),] %>% 
                        st_geometry %>% 
                        st_coordinates %>%
                        st_linestring
    } )
lc_lines = st_sfc(l_lines,
                  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

make_names = function (a, b) {
    paste(
          ifelse(a < b, a, b),
          ifelse(a < b, b, a),
          sep="_"
    )
}

pair_stats = st_sf(
               pairs_data,
               row.names = make_names(
                    stats$short_name[match(pairs_data$loc1, stats$site_name)],
                    stats$short_name[match(pairs_data$loc2, stats$site_name)]
               ),
               geometry = lc_lines
)

pair_stats$distance = st_length(st_geometry(pair_stats))


###
# project points onto the line between the furthest-apart points

north_pt = "conness"
south_pt = "army"
longest = st_geometry(pair_stats)[match("army_conness", rownames(pair_stats))]
flat_coords = 3310
perp_lines = st_nearest_points(st_transform(stats, flat_coords), st_transform(longest, flat_coords))
perp_points = st_cast(perp_lines, "POINT")[2*seq_len(nrow(stats))]
stats$pct = as.vector(st_distance( perp_points, st_transform(stats[south_pt,], flat_coords) ))
stats = stats[order(stats$pct),]

# shows projections -- probably optional
png(file=sprintf("%s.proj.png", basename), width=6.5*288, height=8*288, res=288)
    plot_setup() +
        geom_sf(data=longest, col='blue', lwd=2) +
        geom_sf(data=perp_lines, col='purple') +
        geom_sf(data=perp_points, col='red', pch=20) +
        geom_sf(data=locations, col='black', pch=20)
dev.off()


###
# Heterozygosity

png(file=sprintf("%s.het.png", basename), width=6*288, height=8*288, res=288)

    plot_setup() +
        geom_sf(data=stats, mapping=aes(color=het), cex=6)

dev.off()


###
# IBD

ibdfile = sprintf("%s.ibd.png", basename)

png(file=ibdfile, width=6.5*288, height=8*288, res=288)
layout(1:2)
for (vn in c("dxy", "Fst")) {
    plot(pair_stats[["distance"]], pair_stats[[vn]], pch=20,
         ylab=vn,
         xlab="distance",
         main=basename
    )
}
dev.off()


pdf(file=sprintf("%s.dxy_by_pct.pdf", basename), width=6.5, height=15, pointsize=10)
    yscale = 1.0 * max(pair_stats$dxy, na.rm=TRUE)
    par(mar=c(5, 10, 3, 1)+.1)
    plot(0, type='n',
         xlim=range(stats$pct),
         xlab='PCT distance [m]',
         ylim=c(0, nrow(stats) * yscale), yaxt='n',
         ylab='',
         main="dxy"
    )
    axis(2, las=2, tick=FALSE,
         at=(1:nrow(stats) - 0.5) * yscale,
         labels=rownames(stats)
    )
    abline(h=(0:nrow(stats)) * yscale, lwd=2, col=adjustcolor('black', 0.5))
    for (k in 1:nrow(stats)) {
        ref = rownames(stats)[k]
        ii = match(make_names(ref, rownames(stats)), rownames(pair_stats))
        y = yscale * (k-1) + pair_stats$dxy[ii]
        abline(h=y[k], lty=3, col=adjustcolor(k, 0.5))
        lines(stats$pct, y, col=k, lty=k, type='b', pch=20)
        points(stats$pct[k], y[k], pch="*", cex=5, col=k)
    }
dev.off()

