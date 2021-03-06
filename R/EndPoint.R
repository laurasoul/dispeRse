#' Random walk on a sphere
#'
#' This function allows you calculate the final coordinates and bearing after one step of a random walk on a sphere.
#' @param longitude starting longitude
#' @param latitude starting latitude
#' @param bearing heading in degrees from north
#' @param distance length of step in km
#' @param EarthRad Earth radius in kilometres.
#' @return longitude final longitude
#' @return latitude final latitude
#' @return bearing final bearing in degrees from north
#' @keywords random walk
#' @author Laura C. Soul \email{lauracsoul@@gmail.com}
#' @export
#' @examples
#' EndPoint(longitude = 0, latitude = 0, bearing = 90, distance = 111, EarthRad = 6367.4447)

# Perhaps add conditionals to overcome floating point error when travelling exactly N, E, S, or W
# Similarly perhaps add something to ensure travelling from poles makes sense, e.g. going North from North Pole.

EndPoint <- function(longitude = 0, latitude = 0, bearing = 0, distance = 1, EarthRad = 6367.4447) {
    
    slat_rad <- pi * latitude / 180
    slong_rad <- pi * longitude / 180
    bearing_rad <- pi * bearing / 180

    end_lat_rad <- asin(sin(slat_rad) * cos(distance / EarthRad) + cos(slat_rad) * sin(distance / EarthRad) * cos(bearing_rad))
    end_long_rad <- slong_rad + atan2(sin(bearing_rad) * sin(distance / EarthRad) * cos(slat_rad), cos(distance / EarthRad) - sin(slat_rad) * sin(end_lat_rad))
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

    result <- list(end_long, end_lat, final_bearing)
    names(result) <- c("longitude", "latitude", "bearing")

    return(result)
	
}
