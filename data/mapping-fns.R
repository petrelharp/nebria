library(stars)
library(sf)
library(tidyverse)
library(tigris)

# the PROJ.4 string we're using
.proj4 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

.mapping.thisdir <- file.path(normalizePath("."))

locations <- read.csv(file.path(.mapping.thisdir, "Nebria_ingens_samplesize.csv"), header=TRUE, stringsAsFactors=FALSE) %>%
              st_as_sf(
                       coords = c("longitude", "latitude"), dim="XY",
                       crs = .proj4
                      )

#' Get a slightly expanded SpatialPolygon (rectangle) enclosing the object
#' @param x The spatial object
get_rect <- function (x, fact=.2) {
    the.bbox <- matrix( st_bbox(x)[1:4], nrow=2 )
    the.bbox <- t( apply( the.bbox, 1, function (z) { mean(z) + (1+fact)*(z-mean(z)) } ) )
    st_polygon( list( matrix( c( 
                          the.bbox[c(1,3,3,1,1)],
                          the.bbox[c(2,2,4,4,2)]
                       ), ncol=2 )),
               dim="XY" 
    )
}

#' Get county outlines to overlay on a spatial object
#' @param x The spatial object
get_counties <- function (x, ...) {
    clines = tigris::counties("California", cb=TRUE)
    st_crop(clines, get_rect(x, ...))
}

#' Get state lines to overlay on a spatial object
#' @param x The spatial object
get_statelines <- function (x) {
    tigris::states(c("California", "Oregon"), cb=TRUE)
}

#' Get elevation contours to overlay on a spatial object
#' @param x The spatial object
get_contours <- function (...) {
    st_contour(read_stars("merged_srtm.tif"), breaks=500*(1:8)) # in meters
}

#' Get elevation shading to overlay on a spatial object
#' @param x The spatial object
#' Plot the result with, e.g.
#' plot( shade, col=adjustcolor(grey(seq(0,1,length.out=101)),0.5), legend=FALSE )
get_shading <- function (...) {
    read_stars("cropped_SR_HR.tif")
}


shade = get_shading()
contour = get_contours()

plot_setup <- function (...) {

    ggplot() + 
        geom_stars(data=shade) +
        scale_fill_gradientn(colours=c("#000000BB", "#FFFFFFBB"), guide='none') +
        geom_sf(data=contour, colour=adjustcolor("black", 0.1), fill=NA) # + geom_sf(data=locations)

}

