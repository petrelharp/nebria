#!/usr/bin/env R

args = commandArgs(TRUE)

if (length(args) != 1) {
    stop("Usage: plot_stats.R (name of tree sequence)")
}

basename = gsub(".trees$", "", args[1])

library(png)
library(tidyverse)
library(sf)
library(stars)

source("data/mapping-fns.R", chdir=TRUE)

## too big
# glacier = read_sf("data/glacier_boundary")
suitability = read_stars("geo_only_suitability.tif")

stats_data = read.csv(sprintf("%s.stats.csv", basename))
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

pair_lines = st_sf(
               pairs_data,
               row.names = make_names(
                    stats$short_name[match(pairs_data$loc1, stats$site_name)],
                    stats$short_name[match(pairs_data$loc2, stats$site_name)]
               ),
               geometry = lc_lines
)

pair_lines$distance = st_length(st_geometry(pair_lines))


###
# project points onto the line between the furthest-apart points

north_pt = "conness"
south_pt = "army"
longest = st_geometry(pair_lines)[match("army_conness", rownames(pair_lines))]
flat_coords = 3310
perp_lines = st_nearest_points(st_transform(stats, flat_coords), st_transform(longest, flat_coords))
perp_points = st_cast(perp_lines, "POINT")[2*seq_len(nrow(stats))]
stats$pct = st_distance( perp_points, st_transform(stats[south_pt,], flat_coords) )

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
        geom_sf(data=stats, mapping=aes(color=sim_pi), cex=6)

dev.off()


###
# IBD

ibdfile = sprintf("%s.ibd.png", basename)

png(file=ibdfile, width=6.5*288, height=8*288, res=288)
layout(1:2)
for (vn in c("dxy", "Fst")) {
    plot(pair_lines[["distance"]], pair_lines[[vn]], pch=20,
         ylab=vn,
         xlab="distance",
         main=basename
    )
}
dev.off()


k = 1
ref = stats[k,]


plot_setup() +
    geom_sf(data=stats, mapping=aes(col=pct), cex=5, pch=20)

