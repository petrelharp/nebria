library(raster)
library(png)

# proportion of pixels which should not be zeroed out
keep_prop <- 1/25

set.seed(123)

fn <- "geo_only_suitability"
x <- raster(paste0(fn, ".tif"))
png(file=paste0(fn, "_with_axes.png"),
    width=48*dim(x)[1]*6/max(dim(x)),
    height=48*dim(x)[2]*6/max(dim(x)),
    res=48,
    pointsize=10)
raster::plot(x)
dev.off()
png::writePNG(as.matrix(x), paste0(fn, ".png"), dpi=24)

# Project all to align with "current", which is a factor of 5
tifs <- c("current.tif", "LH0_4.tif", "MH4_8.tif", "EH8_12.tif", "BA13_15.tif", "HS15_17.tif")
rasters <- lapply(tifs, raster)
names(rasters) <- gsub(".tif$", "", tifs)
for (k in 2:length(tifs)) {
    rasters[[k]] <- projectRaster(rasters[[k]], rasters[['current']])
    stopifnot(res(rasters[[k]]) == res(rasters[['current']]))
}

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


### choose actually good patches, which are common across time periods
keep <- rasters[['current']]
values(keep) <- ifelse(runif(prod(dim(keep))) < keep_prop, 1.0, 0.0)

for (k in seq_along(rasters)) {
    fn <- names(rasters)[k]
    x <- rasters[[k]]
    # version with axes
    png(file=paste0(fn, "_with_axes.png"),
        width=48*dim(x)[1]*6/max(dim(x)),
        height=48*dim(x)[2]*6/max(dim(x)),
        res=48,
        pointsize=10)
    raster::plot(x)
    dev.off()

    # plain png version: 
    # R = all (original) habitat
    # G = subset of habitat
    # B = unused
    x[is.na(x)] <- 0.0
    xm <- as.matrix(x)
    xmk <- as.matrix(x * keep)
    am <- array(c(xm, xmk, xm * 0.0), dim=c(dim(xm), 3))
    png::writePNG(am, paste0(fn, ".png"), dpi=24)
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

fn <- "geo_only_suitability"
x <- crop(raster(paste0(fn, ".tif")), small_box)
writeRaster(x, file.path(outdir, paste0(fn, ".tif")))
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

