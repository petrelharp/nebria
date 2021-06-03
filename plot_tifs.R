library(raster)

tifs <- list.files(".", "*.tif")

for (fn in tifs) {
    x <- raster(fn)
    png(file=gsub("tif", "png", fn),
        width=48*dim(x)[1]*6/max(dim(x)),
        height=48*dim(x)[2]*6/max(dim(x)),
        res=48,
        pointsize=10)
    plot(x)
    dev.off()
}

