#' Random walk on a sphere
#'
#' This function allows you calculate the final coordinates and bearing after one step of a random walk on a sphere.
#' @param latitude starting latitude
#' @param longitude starting longitude
#' @param bearing heading in degrees from north
#' @param distance length of step in km
#' @return lat final latitude
#' @return long final longitude
#' @return bearing final bearing in degrees from north
#' @keywords random walk
#' @export
#' @examples
#' endpoint(slat=0, slong=0, bearing=90, distance=111)

# Perhaps add conditionals to overcome floating point error when travelling exactly N, E, S, or W
# Similarly perhaps add something to ensure travelling from poles makes sense, e.g. going North from North Pole.

endpoint <- function(slat=0, slong=0, bearing=0, distance=1) {
    
    R <- 6367.4447 #average radius

    slat_rad <- pi * slat / 180
    slong_rad <- pi * slong / 180
    bearing_rad <- pi * bearing / 180

    end_lat_rad <- asin(sin(slat_rad) * cos(distance / R) + cos(slat_rad) * sin(distance / R) * cos(bearing_rad))
    end_long_rad <- slong_rad + atan2(sin(bearing_rad) * sin(distance / R) * cos(slat_rad), cos(distance / R) - sin(slat_rad) * sin(end_lat_rad))
    end_long_rad <- (end_long_rad + pi) %% (2 * pi) - pi

    end_lat <- end_lat_rad * 180 / pi
    end_long <- end_long_rad * 180 / pi

    dLon <- slong_rad - end_long_rad
    y <- sin(dLon) * cos(slat_rad)
    x <- cos(end_lat_rad) * sin(slat_rad) - sin(end_lat_rad) * cos(slat_rad) * cos(dLon)
    z <- atan2(y, x)
    final_bearing_rad <- z + pi
    final_bearing <- final_bearing_rad * 180 / pi
    final_bearing <- (final_bearing + 360) %% 360

    result <- list(end_lat, end_long, final_bearing)
    names(result) <- c("lat", "long", "bearing")

    return(result)
	
}
