source("mapping-fns.R", chdir=TRUE)
library(raster)

e = raster::extent(as.vector(st_bbox(get_rect(locations, 0.8))[c("xmin", "xmax", "ymin", "ymax")]))

# DEM

dems = list.files("data", "^srtm_.*tif", full.names=TRUE)

# from https://stackoverflow.com/questions/15876591/merging-multiple-rasters-in-r

template = raster(e, resolution=c( 0.0008333333, 0.0008333333))
raster::projection(template) <- .proj4

dem = crop(raster(dems[1]), template)
template = projectRaster(from=dem, to=template, alignOnly=TRUE)
for (demf in dems) {
    dem = raster(demf)
    template = mosaic(template, crop(dem, template), fun=mean)
}
template = raster::aggregate(template, fact=5)

raster::writeRaster(template, file="merged_srtm.tif", format="GTiff", overwrite=TRUE)


# Shading

template = raster(e, resolution=c( 0.0008333333, 0.0008333333))
raster::projection(template) <- .proj4
msr <- projectRaster( raster("US_MSR_10M/US_MSR.tif"), crs=crs(template))
writeRaster( crop( msr, template ), file="cropped_SR_HR.tif", overwrite=TRUE)


