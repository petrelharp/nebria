library(raster)
library(sp)
library(sf)
library(png)

# proportion of pixels which should not be zeroed out
keep_prop <- 1/8
suitability <- "Geo_only_slope_aspect_drainage"

set.seed(123)

sample_locs <- read.csv("data/Nebria_ingens_samplesize.csv")
samples <- SpatialPointsDataFrame(
                  coords = sample_locs[c("longitude", "latitude")],
                  data = sample_locs[setdiff(names(sample_locs), c("longitude", "latitude"))],
                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
)


for (fn in file.path("geo_layers", c(suitability, "glacier_boundary_fill"))) {
    x <- raster(paste0(fn, ".tif"))
    png(file=paste0(fn, "_with_axes.png"),
        width=48*dim(x)[1]*6/max(dim(x)),
        height=48*dim(x)[2]*6/max(dim(x)),
        res=48,
        pointsize=10)
    raster::plot(x)
    points(samples)
    dev.off()
    png::writePNG(as.matrix(x), paste0(fn, ".png"), dpi=24)
}

# Project all to align with "current", which is a factor of 5
tifs <- file.path("geo_layers",
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
names(rasters) <- basename(gsub(".tif$", "", tifs))
for (k in 2:length(tifs)) {
    rasters[[k]] <- projectRaster(rasters[[k]], rasters[['current']])
    stopifnot(res(rasters[[k]]) == res(rasters[['current']]))
}

# Compute spatial resolution of grid cells in km,
# to input into the MAP_RES parameter in SLiM
xy <- SpatialPoints(cbind(
                c(-119.5, -118.5),
                c(37.0, 37.0)),
            proj4string = CRS("+proj=longlat"))
xy_coords <- xyFromCell(rasters[["current"]], cellFromXY(rasters[["current"]], xy))
ij <- rowColFromCell(rasters[["current"]], cellFromXY(rasters[["current"]], xy_coords))
xy <- SpatialPoints(xy_coords, proj4string = CRS(proj4string(rasters[["current"]])))
dist_xy <- pointDistance(xy)[2,1]  # in meters
dist_ij <- sqrt((ij[1,1] - ij[2,1])^2 + (ij[1,2] - ij[2,2])^2)
res_horiz <- dist_xy / dist_ij

# and now vertically
xy <- SpatialPoints(cbind(
                c(-119.0, -119.0),
                c(36.1, 37.9)),
            proj4string = CRS("+proj=longlat"))
xy_coords <- xyFromCell(rasters[["current"]], cellFromXY(rasters[["current"]], xy))
ij <- rowColFromCell(rasters[["current"]], cellFromXY(rasters[["current"]], xy_coords))
xy <- SpatialPoints(xy_coords, proj4string = CRS(proj4string(rasters[["current"]])))
dist_xy <- pointDistance(xy)[2,1]  # in meters
dist_ij <- sqrt((ij[1,1] - ij[2,1])^2 + (ij[1,2] - ij[2,2])^2)
res_vert <- dist_xy / dist_ij

cat(sprintf("The resolutions of the 'current' raster are %f (E-W) and %f (N-S) meters per pixel.\n", res_horiz, res_vert))


# Convert sample locations to coordinates in SLiM,
# which will be in units of kilometers from the *center* of the pixel in the SW corner of the map
xy_to_slim <- function (xy) {
    # verified that xyFromCell finds the *center* of the cell
    ll_xy <- c(extent(rasters[['current']])[1],
               extent(rasters[['current']])[3])
    ll_cell <- cellFromXY(rasters[['current']], ll_xy)
    ll_xy <- SpatialPoints(
                    xyFromCell(rasters[['current']], ll_cell),
                    proj4string=CRS(proj4string(rasters[['current']]))
            )
    triangle_point <- SpatialPoints(
                    cbind(
                        rep(coordinates(ll_xy)[1], length(xy)),
                        coordinates(xy)[,2]
                    ))
    horiz_dist <- pointDistance(triangle_point, xy, lonlat=TRUE) / 1000
    vert_dist <- pointDistance(ll_xy, triangle_point, lonlat=TRUE) / 1000
    return(data.frame(slim_x=horiz_dist, slim_y=vert_dist))
}
sample_locs <- cbind(sample_locs, xy_to_slim(samples))
write.csv(sample_locs, "sample_locs.csv", row.names=FALSE)

# Fit a linear model using xy_to_slim to convert slim coordinates back to latitude longitude
box <- c(-119.3971, -118.2471, 36.41429, 37.94329)
input_xy <- expand.grid(long = seq(from = box[1], to = box[2], length.out = 10), lat = seq(from = box[3], to = box[4], length.out = 10))
input_xy <- cbind(input_xy, xy_to_slim(input_xy))
lat_model <- lm(lat ~ slim_y, data = input_xy)
long_model <- lm(long ~ slim_x*slim_y, data = input_xy)

# Check model

#test <- cbind(input_xy, slim_to_latlong(input_xy))

# cut out the white mountains
wm_box <- SpatialPolygons(list(Polygons(list(
              Polygon(rbind(c(-118, 37.3),
                            c(-118, 38),
                            c(-118.5, 38),
                            c(-118.5, 37.3)))),
            ID="wm")),
            proj4string=CRS(proj4string(rasters[['current']]))
)
rasters <- lapply(rasters, mask, wm_box, inverse=TRUE)

# remove the glacier from the oldest time
glacier_mask <- 
    projectRaster(
        raster("geo_layers/glacier_boundary_fill.tif"),
        rasters[["current"]],
        method='ngb',
    )
rasters[["21000_LGM_CCSM"]] <- mask(rasters[["21000_LGM_CCSM"]], glacier_mask, maskvalue=2, updatevalue=0.0)

### choose actually good patches, which are common across time periods
keep <- rasters[['current']]
values(keep) <- ifelse(runif(prod(dim(keep))) < keep_prop, 1.0, 0.0)

for (k in seq_along(rasters)) {
    fn <- names(rasters)[k]
    x <- rasters[[k]]
    # version with axes
    png(file=file.path("geo_layers", paste0(fn, "_with_axes.png")),
        width=144*dim(x)[1]*12/max(dim(x)),
        height=144*dim(x)[2]*12/max(dim(x)),
        res=48,
        pointsize=10)
    raster::plot(x)
    points(samples)
    dev.off()

    # plain png version: 
    # R = all (original) habitat
    # G = subset of habitat
    # B = unused
    x[is.na(x)] <- 0.0
    xm <- as.matrix(x)
    xmk <- as.matrix(x * keep)
    am <- array(c(xm, xmk, xm * 0.0), dim=c(dim(xm), 3))
    png::writePNG(am, file.path("geo_layers", paste0(fn, ".png")), dpi=24)
}


### also write out a smaller example
outdir <- "small"
if (!file.exists(outdir)) {
    dir.create(outdir)
}
small_box <- SpatialPolygons(list(Polygons(list(
              Polygon(rbind(c(-118.4, 37.2),
                            c(-118.4, 37),
                            c(-118.6, 37),
                            c(-118.6, 37.2)))),
            ID="small")),
            proj4string=CRS(proj4string(rasters[['current']]))
)

fn <- suitability
x <- crop(raster(file.path("geo_layers", paste0(fn, ".tif"))), small_box)
writeRaster(x, file.path(outdir, paste0(fn, ".tif")), overwrite=TRUE)
png::writePNG(as.matrix(x), file.path(outdir, paste0(fn, ".png")), dpi=24)

for (k in seq_along(rasters)) {
    fn <- names(rasters)[k]
    x <- rasters[[k]]
    x[is.na(x)] <- 0.0
    xm <- as.matrix(crop(x, small_box))
    xmk <- as.matrix(crop(x * keep, small_box))
    am <- array(c(xm, xmk, xm * 0.0), dim=c(dim(xm), 3))
    png::writePNG(am, file.path(outdir, paste0(fn, ".png")), dpi=24)
}



### and, yosemite
outdir <- "yosemite"
if (!file.exists(outdir)) {
    dir.create(outdir)
}
yo_box <- SpatialPolygons(list(Polygons(list(
              Polygon(rbind(c(-119.07, 38),
                            c(-119.07, 37.5),
                            c(-119.6, 37.5),
                            c(-119.6, 38)))),
            ID="yosemite")),
            proj4string=CRS(proj4string(rasters[['current']]))
)

fn <- suitability
x <- crop(raster(file.path("geo_layers", paste0(fn, ".tif"))), small_box)
writeRaster(x, file.path(outdir, paste0(fn, ".tif")), overwrite=TRUE)
png::writePNG(as.matrix(x), file.path(outdir, paste0(fn, ".png")), dpi=24)

for (k in seq_along(rasters)) {
    fn <- names(rasters)[k]
    x <- rasters[[k]]
    x[is.na(x)] <- 0.0
    xm <- as.matrix(crop(x, small_box))
    xmk <- as.matrix(crop(x * keep, small_box))
    am <- array(c(xm, xmk, xm * 0.0), dim=c(dim(xm), 3))
    png::writePNG(am, file.path(outdir, paste0(fn, ".png")), dpi=24)
}


