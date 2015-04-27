#' One-to-many Great Circle Distances
#'
#' Returns a set of distances from a single point of longitude and latitude
#' @param point_longitude Decimalised longitude of single point.
#' @param point_latitude Decimalised latitude of single point.
#' @param point_longitudes Decimalised longitudes of many points.
#' @param point_latitudes Decimalised latitudes of many points.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return A vector of minimum Great Circle Distances between single point and each of many points.
#' @details Nothing yet.
#' @examples
#' One2ManyGreatCircleDistance(0, 0, c(-179, 0, 179), c(89, 0, -89))

One2ManyGreatCircleDistance <- function(point_longitude, point_latitude, point_longitudes, point_latitudes, EarthRad = 6367.4447) {
	
	# Create empty square matrix:
	GC_distances <- vector(mode="numeric")
	
	# For each first coordinate in pariwise comparison:
	for(i in 1:length(point_longitudes)) {
		
		# Store minimal Great Circle distance:
		GC_distances[i] <- GreatCircleDistanceFromLongLat(point_longitude, point_latitude, point_longitudes[i], point_latitudes[i], EarthRad = EarthRad)
		
	}
	
	# Return Great Circle distances:
	return(GC_distances)
	
}
