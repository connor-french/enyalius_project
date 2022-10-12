library(raster)
library(sp)
library(sf)
library(png)
library(here)


# Compute spatial resolution of grid cells in km,
# to input into the MAP_RES parameter in SLiM
get_slim_res <- function(rast) {
  e <- extent(rast)

  xy <- SpatialPoints(cbind(
    c(e@xmin, e@xmax),
    c(e@ymin, e@ymin)),
    proj4string = CRS("+proj=longlat"))

  xy_coords <- xyFromCell(rast, cellFromXY(rast, xy))
  ij <- rowColFromCell(rast, cellFromXY(rast, xy_coords))
  xy <- SpatialPoints(xy_coords, proj4string = CRS(proj4string(rast)))
  dist_xy <- pointDistance(xy)[2,1]  # in meters
  dist_ij <- sqrt((ij[1,1] - ij[2,1])^2 + (ij[1,2] - ij[2,2])^2)
  res_horiz <- dist_xy / dist_ij

  xy <- SpatialPoints(cbind(
    c(e@xmin, e@xmin),
    c(e@ymin, e@ymax)),
    proj4string = CRS("+proj=longlat"))

  xy_coords <- xyFromCell(rast, cellFromXY(rast, xy))
  ij <- rowColFromCell(rast, cellFromXY(rast, xy_coords))
  xy <- SpatialPoints(xy_coords, proj4string = CRS(proj4string(rast)))
  dist_xy <- pointDistance(xy)[2,1]  # in meters
  dist_ij <- sqrt((ij[1,1] - ij[2,1])^2 + (ij[1,2] - ij[2,2])^2)
  res_vert <- dist_xy / dist_ij

  return(c(res_horiz, res_vert))
}



# Convert sample locations to coordinates in SLiM,
# which will be in units of kilometers from the *center* of the pixel in the SW corner of the map
# assumes coordinates are not projected

xy_to_slim <- function (locs, rast) {

  xy <- SpatialPointsDataFrame(
    coords = locs[c("longitude", "latitude")],
    data = locs[setdiff(names(locs), c("longitude", "latitude"))],
    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  )

  # verified that xyFromCell finds the *center* of the cell
  ll_xy <- c(extent(rast)[1],
             extent(rast)[3])
  ll_cell <- cellFromXY(rast, ll_xy)
  ll_xy <- SpatialPoints(
    xyFromCell(rast, ll_cell),
    proj4string=CRS(proj4string(rast))
  )
  triangle_point <- SpatialPoints(
    cbind(
      rep(coordinates(ll_xy)[1], length(xy)),
      coordinates(xy)[,2]
    ))
  horiz_dist <- pointDistance(triangle_point, xy, lonlat=TRUE) / 1000

  vert_dist <- pointDistance(ll_xy, triangle_point, lonlat=TRUE) / 1000

  sample_locs <- cbind(locs, data.frame(slim_x=horiz_dist, slim_y=vert_dist))

  return(sample_locs)
}


#### Not so sure what the purpose of this is. I'm going to do it my way first with "write_patches_slim"
#### Keep cells that are consistently occupied across time periods

keep_patches_slim <- function(keep_prop, current_rast, rast_stack, locs) {
  keep <- current_rast
  values(keep) <- ifelse(runif(prod(dim(keep))) < keep_prop, 1.0, 0.0)

  xy <- SpatialPointsDataFrame(
    coords = locs[c("longitude", "latitude")],
    data = locs[setdiff(names(locs), c("longitude", "latitude"))],
    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  )

  for (k in seq_along(rast_stack)) {
    fn <- names(rast_stack)[k]
    x <- rast_stack[[k]]
    # version with axes
    png(file= here("analysis", "slim", "geo_layers", paste0(fn, "_with_axes.png")),
        width=144*dim(x)[1]*12/max(dim(x)),
        height=144*dim(x)[2]*12/max(dim(x)),
        res=48,
        pointsize=10)
    raster::plot(x)
    points(xy)
    dev.off()

    # plain png version:
    # R = all (original) habitat
    # G = subset of habitat
    # B = unused
    x[is.na(x)] <- 0.0
    xm <- as.matrix(x)
    xmk <- as.matrix(x * keep)
    am <- array(c(xm, xmk, xm * 0.0), dim=c(dim(xm), 3))
    png::writePNG(am, here("analysis", "slim", "geo_layers", paste0(fn, "_with_axes.png")), dpi=24)
  }

}


write_patches_slim <- function(current_rast, rast_stack=NULL, locs) {


  xy <- SpatialPointsDataFrame(
    coords = locs[c("longitude", "latitude")],
    data = locs[setdiff(names(locs), c("longitude", "latitude"))],
    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  )

  if (is.null(rast_stack)) {
    # plain png version:
    # R = all (original) habitat
    # G = subset of habitat
    # B = unused
    current_rast[is.na(current_rast)] <- 0.0
    xm <- as.matrix(current_rast)
    xmk <- as.matrix(current_rast * current_rast)
    am <- array(c(xm, xmk, xm * 0.0), dim=c(dim(xm), 3))
    png::writePNG(am, here("analysis", "slim", "slim_layers", paste0(names(current_rast), ".png")), dpi=24)
  } else {
    for (k in seq_along(rast_stack)) {
      fn <- names(rast_stack)[k]
      x <- rast_stack[[k]]
      # version with axes
      png(file= here("analysis", "slim", "slim_layers", paste0(fn, "_with_axes.png")),
          width=144*dim(x)[1]*12/max(dim(x)),
          height=144*dim(x)[2]*12/max(dim(x)),
          res=48,
          pointsize=10)
      raster::plot(x)
      points(xy)
      dev.off()

      # plain png version:
      # R = all (original) habitat
      # G = subset of habitat
      # B = unused
      x[is.na(x)] <- 0.0
      xm <- as.matrix(x)
      xmk <- as.matrix(x * current_rast)
      am <- array(c(xm, xmk, xm * 0.0), dim=c(dim(xm), 3))
      png::writePNG(am, here("analysis", "slim", "slim_layers", paste0(fn, ".png")), dpi=24)
    }
  }



}

