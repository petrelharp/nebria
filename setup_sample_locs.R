#!/usr/bin/env R

# This just adds the PCT coordinate to the "sample_locs.csv" file.

library(sf)
library(stars)

sample_locs = read.csv("sample_locs.csv")

locs = st_as_sf(
            sample_locs, agr="identity",
            coords = c("longitude", "latitude"), dim="XY",
            crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
)

north_pt = "Conness Lake"
south_pt = "Army Pass"

longest <- st_sfc(
              locs[c(match(south_pt, locs$site_name), 
                  match(north_pt, locs$site_name)),] %>%
               st_geometry %>%
               st_coordinates %>%
               st_linestring,
          crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

flat_coords = 3310
perp_lines = st_nearest_points(st_transform(locs, flat_coords), st_transform(longest, flat_coords))
perp_points = st_cast(perp_lines, "POINT")[2*seq_len(nrow(locs))]
sample_locs$pct = as.vector(st_distance( perp_points, st_transform(locs[match(south_pt, locs$site_name),], flat_coords) ))

write.csv(sample_locs, file="sample_locs.csv", row.names=FALSE)
