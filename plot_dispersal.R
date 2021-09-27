library(raster)

suitability <- raster::raster("small/geo_only_suitability.tif")
habitat <- png::readPNG("small/BA13_15.png")
all_habitat <- t(habitat[rev(seq_len(nrow(habitat))),,1])
habitat <- t(habitat[rev(seq_len(nrow(habitat))),,2])

defineConstant <- assign
defineConstant("P_D", 0.1); # // dispersal probability
defineConstant("MEAN_HOPS", 20); # // mean number of dispersal hops
defineConstant("DISPERSAL_SIGMA", 1.0); # // dispersal distance
defineConstant("MAP_RES_X", 0.742); # // resolution in km (see make_maps.R)
defineConstant("MAP_RES_Y", 0.925); # // resolution in km (see make_maps.R)
defineConstant("STEPSIZE", 1); #   // years per generation (to make debugging go faster)
defineConstant("SIGMA", DISPERSAL_SIGMA * sqrt(STEPSIZE));

stopifnot(STEPSIZE == 1)

loc_to_map <- function (loc) {
    loc <- rbind(loc)
    x0 <- extent(suitability)[1]
    x1 <- extent(suitability)[2]
    y0 <- extent(suitability)[3]
    y1 <- extent(suitability)[4]
    x <- x0 + (x1 - x0) * loc[,1] / MAP_RES_X / nrow(habitat)
    y <- y0 + (y1 - y0) * loc[,2] / MAP_RES_Y / ncol(habitat)
    return(cbind(x,y))
}

p1.pointInBounds <- function (loc) {
    xy <- loc_to_map(loc)
    bounds <- extent(suitability)
    return(
           xy[,1] >= bounds[1]
           & xy[,1] <= bounds[2]
           & xy[,2] >= bounds[3]
           & xy[,2] <= bounds[4]
    )
}

p1.spatialMapValue <- function(map, loc) {
    if (p1.pointInBounds(loc)) {
        xy <- loc_to_map(loc)
        out <- raster::extract(get(map), xy)
    } else {
        out <- 0.0
    }
    return(out)
}

reproduction <- function (loc, P_D=1) {
    path <- rbind(c(loc, NA))
            loc = loc + rnorm(2, mean=0, sd=SIGMA / sqrt(MEAN_HOPS));
    path <- rbind(c(loc, 1))
			num_hops = 1;
			tries = 0;
			if (runif(1) < P_D) {
                # // NegBinom(r, p) is Pois(Gamma(shape=r, rate=(1-p)/p))
                k = 1 + rpois(1, rexp(1, rate=1/MEAN_HOPS));
				while (num_hops < k & tries < 100) {
                    current = p1.spatialMapValue('suitability', loc);
					next_loc = loc + rnorm(2, mean=0, sd=SIGMA / sqrt(MEAN_HOPS));
					if (runif(1) * current <= p1.spatialMapValue('suitability', next_loc)) {
						loc = next_loc;
						num_hops = num_hops + 1;
    path <- rbind(path, c(next_loc, 1)) } else { path <- rbind(path, c(next_loc, 0))
					}
					tries = tries + 1;
				}
			}
    stopifnot(sum(path[,3], na.rm=TRUE) == num_hops)
    stopifnot(nrow(path) == tries + 1)
    return(path)
}

plot_path <- function (path) {
    k0 <- rep(NA, nrow(path) - 1)
    k1 <- rep(NA, nrow(path) - 1)
    last_k <- 1
    for (k in seq_len(NROW(path))[-1]) {
        k0[k] <- last_k
        k1[k] <- k
        if (path[k, 3] > 0) {
            last_k <- k
        }
    }
    xy0 <- loc_to_map(path[k0, 1:2])
    xy1 <- loc_to_map(path[k1, 1:2])
    segments(
        x0=xy0[,1],
        y0=xy0[,2],
        x1=xy1[,1],
        y1=xy1[,2],
        col=ifelse(path[,3] == 0, "red", "black")
    )
}

good_locs <- cbind(row(habitat)[habitat > 0.1] * MAP_RES_X,
                   col(habitat)[habitat > 0.1] * MAP_RES_Y)
stopifnot(all(p1.pointInBounds(good_locs)))

good_points <- loc_to_map(good_locs)

paths <- lapply(1:nrow(good_locs), function (k)
                reproduction(good_locs[k,])
            )

# note NegBinom(r, p) is Pois(Gamma(shape=r, rate=(1-p)/p))
rdispersal_dist <- function (n, shape=100) {
    nhops <- 1 + rpois(n, rgamma(n, shape=shape, rate=shape/MEAN_HOPS))
    xy <- cbind(
            rnorm(n, sd=sqrt(nhops) * SIGMA / sqrt(MEAN_HOPS)),
            rnorm(n, sd=sqrt(nhops) * SIGMA / sqrt(MEAN_HOPS))
    )
    return(sqrt(rowSums(xy^2)))
}

pdf(file="example_dispersal.pdf", width=10, height=5, pointsize=10)
layout(t(1:3))
    x <- rdispersal_dist(1e4, shape=1)
    hist(x, breaks=100,
         main="dispersal distance (km)")
    xx <- seq(0.1, 10, length.out=101)
    plot(xx, type='l',
         sapply(xx, function (z) mean(x * z / SIGMA > 5)),
         xlab='SIGMA',
         ylab='proportion > 5',
         log='y'
    )
    abline(v=1, lty=3)
    plot(suitability, main=sprintf("SIGMA=%0.1f", SIGMA))
    points(good_points, cex=5)
    for (path in paths) {
        plot_path(path)
    }
    scalebar(type='bar', lonlat=TRUE, below='kilometers')
dev.off()

