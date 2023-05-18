library(tidyverse)
library(stars)
library(sf)
library(tigris)
library(raster)
library(ggforce)
library(viridis)

# Convert SLiM locations to coordinates
raster_bbox <- c(-119.9996, -117.99993, 35.9995, 38.01215)
res_x <- 0.742 # resolution in km (see make_maps.R)
res_y <- 0.925 #resolution in km (see make_maps.R)
horizontal_distance <- pointDistance(raster_bbox[c(1, 3)], raster_bbox[c(2, 3)], lonlat = TRUE)
vertical_distance <- pointDistance(raster_bbox[c(1, 3)], raster_bbox[c(1, 4)], lonlat = TRUE)
h_slim <- horizontal_distance/742
v_slim <- vertical_distance/925

anc <- read_csv("data_for_ancestry_plot.csv")
all_inds <-  read_csv("data_for_all_inds.csv")
pop_size <- read_csv("../param_grid/post_21000_2022-12-05/run_2022-12-05_000018/sim_1876298096991.log",
                       skip = 1)
ggplot(pop_size, aes(x = years_ago, y = num_individuals)) +
  geom_point() +
  scale_x_reverse() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  xlab("Years ago") +
  ylab("Population size")
ggsave("popsize.png")
ggplot(pop_size, aes(x = years_ago, y = num_patches)) +
  geom_point() +
  scale_x_reverse() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  xlab("Years ago") +
  ylab("Occupied patches")
ggsave("patches.png")
#Remove points with 0 ancestry proportions
anc <- filter(anc, anc_prop >0)
# Add names of individuals
anc <- mutate(anc, ind_name = case_when(ind == 2166 ~ "Northernmost individual", ind == 4299 ~ "Southernmost individual"))
# Add names for time ago
anc <- mutate(anc, time_name = factor(paste0(time, " ya"), levels = c("0 ya", "2250 ya", "6263 ya", "10013 ya",
                                          "12300 ya", "13800 ya", "15850 ya", "20999 ya")))
anc_south <- filter(anc, ind_name == "Southernmost individual")
anc_north <- filter(anc, ind_name == "Northernmost individual")

ggplot(anc_south, aes(x = loc1, y = loc2, size = anc_prop)) + 
  geom_point() +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.title=element_blank()) +
  coord_fixed(ratio = res_y/res_x) +
  xlim(0, h_slim) +
  ylim(0, v_slim) +
  facet_wrap(~time_name) +
  theme(legend.position = "none", strip.text.x = element_text(size = 15))
ggsave("anc_south.png")

ggplot(anc_north, aes(x = loc1, y = loc2, size = anc_prop)) + 
  geom_point() +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.title=element_blank()) +
  coord_fixed(ratio = res_y/res_x) +
  xlim(0, h_slim) +
  ylim(0, v_slim) +
  facet_wrap(~time_name) +
  theme(legend.position = "none", strip.text.x = element_text(size = 15))
ggsave("anc_north.png")

all_inds <- mutate(all_inds, time_name = factor(paste0(time, " ya"), 
                                                levels = c("20999 ya", "15850 ya",
                                                           "13800 ya", "12300 ya",
                                                           "10013 ya", "6263 ya",
                                                           "2250 ya", "0 ya")))
for(j in 1:length(levels(all_inds$time_name))){
  print(j)
  p <- ggplot(all_inds, aes(x = loc1, y = loc2)) +
    geom_point(size =0.1) +
    theme_bw() +
    theme(axis.text=element_blank(),
          axis.title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_fixed(ratio = res_y/res_x) +
    xlim(0, h_slim) +
    ylim(0, v_slim) +
    facet_wrap_paginate(~time_name, nrow = 1, ncol = 1, page = j) +
    theme(strip.text.x = element_text(size = 25),
          strip.background.x = element_blank())
  ggsave(paste0("locations", j, ".png"), plot = p)
}

ggplot(all_inds, aes(x = loc1, y = loc2)) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(axis.text=element_blank(),
          axis.title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_fixed(ratio = res_y/res_x) +
    xlim(0, h_slim) +
    ylim(0, v_slim) +
    facet_wrap(~time_name, nrow = 2) +
    theme(strip.text.x = element_text(size = 15),
          strip.background.x = element_blank())
ggsave("locations.png")


locations <- read.csv("../data/Nebria_ingens_samplesize.csv", header=TRUE, stringsAsFactors=FALSE) %>%
  st_as_sf(
  coords = c("longitude", "latitude"), dim="XY",
  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
)


contour <- st_contour(read_stars("../data/merged_srtm.tif"), breaks=500*(1:8))
shade <- read_stars("../data/cropped_SR_HR.tif")

ggplot() +
  geom_stars(data=shade) +
  scale_fill_gradientn(colours=c("#000000BB", "#FFFFFFBB"), guide='none') +
  geom_sf(data=contour, colour=adjustcolor("black", 0.1), fill=NA)
  geom_sf(data=locations)
ggsave("contours.png")

library("rnaturalearth")
library("rnaturalearthdata")
library(maps)
# Define box that is the sampling area
contour_bbox <- st_bbox(contour)
box_mat <- matrix(contour_bbox[c(1, 2, 1, 4, 3, 4, 3, 2, 1, 2)], ncol = 2, byrow = TRUE)
sample_rect <- st_sfc(st_polygon(list(box_mat)), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

ggplot(data = states) +
  geom_sf() +
  geom_sf(data = sample_rect, fill = NA,linewidth = 0.5) +
  geom_sf(data = locations) +
  coord_sf(xlim = c(-124.5000, -114.99993), ylim = c(32.9995, 45.01215), expand = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("sample_area_big.png", width = 10)

ggplot(data = states) +
  geom_sf() +
  geom_sf(data = sample_rect, fill = NA, linewidth = 0.5) +
  coord_sf(xlim = c(-124.5000, -114.99993), ylim = c(32.9995, 45.01215), expand = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("sample_area_no_locations.png", width = 3, height = 4)

#Species distribution models
tifs <- file.path("../geo_layers",
                  c("current.tif",
                    "00300-04200_LH.tif",
                    "04200_08326_MH.tif",
                    "08326-11700_EH.tif",
                    "11700_12900_YDS.tif",
                    "12900_14700_BA.tif",
                    "14700_17000_HS.tif",
                    "21000_LGM_CCSM.tif")
)
rasters <- lapply(tifs, raster)
rasters_pts <- lapply(rasters, rasterToPoints, spatial = TRUE, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
# Then to a 'conventional' dataframe
rasters_df  <- lapply(rasters_pts, data.frame)
names(rasters_df) <- c("0 ya",
                       "00300-04200 ya",
                       "04200-08326 ya",
                       "08326-11700 ya",
                       "11700-12900 ya",
                       "12900-14700 ya",
                       "14700-17000 ya",
                       "21000 ya")
rename_raster_col <- function(df){
  names(df) <- c("quality", "x", "y", "optional")
  return(df)
}
rasters_df <- lapply(rasters_df, rename_raster_col)
all_rasters <- bind_rows(rasters_df, .id = "time")
all_rasters$time <- factor(all_rasters$time, levels = c("21000 ya",
                                                        "14700-17000 ya",
                                                        "12900-14700 ya",
                                                        "11700-12900 ya",
                                                        "08326-11700 ya",
                                                        "04200-08326 ya",
                                                        "00300-04200 ya",
                                                        "0 ya"))
ggplot() +
  geom_raster(data = all_rasters, aes(x = x, y = y, fill = quality)) +
  geom_sf() +
  facet_wrap(~time, nrow = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(name = "Habitat quality", low = "white", high = "red")
  scale_fill_viridis(name = "Habitat quality", option = "rocket", direction = -1, )
ggsave("sdm.png", width = 10)

test <- data.frame(x = 1:16, y = 1:16, quality = 0:15)
ggplot() +
  geom_raster(data = test, aes(x = x, y = y, fill = quality)) +
  geom_sf() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_viridis(name = "Habitat quality", option = "rocket", direction = -1, begin = )
  scale_fill_gradient(name = "Habitat quality", low = "white", high = "red")
