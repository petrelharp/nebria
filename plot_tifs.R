library(raster)
library(png)

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

# proportion of pixels which should not be zeroed out
keep_prop <- 1/25

# Project all to align with "current", which is a factor of 5
tifs <- c("current.tif", "LH0_4.tif", "MH4_8.tif", "EH8_12.tif", "BA13_15.tif", "HS15_17.tif")
rasters <- lapply(tifs, raster)
names(rasters) <- gsub(".tif$", "", tifs)
for (k in 2:length(tifs)) {
    rasters[[k]] <- projectRaster(rasters[[k]], rasters[['current']])
    stopifnot(res(rasters[[k]]) == res(rasters[['current']]))
}

### choose actually good patches, which are common across time periods
keep <- as.matrix(rasters[['current']])
keep[] <- ifelse(runif(prod(dim(keep))) < keep_prop, 1.0, 0.0)

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
    xm <- as.matrix(x)
    xm[is.na(xm)] <- 0.0
    am <- array(c(xm, xm * keep, xm * 0.0), dim=c(dim(xm), 3))
    png::writePNG(am, paste0(fn, ".png"), dpi=24)
}

