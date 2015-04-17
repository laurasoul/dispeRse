#' Returns great circle distance from a set of coordinates
#'
#' Returns great circle distance in kilometres from a two longitude-latitude points
#'
#' @param long1 Decimalised longitude for first point.
#' @param lat1 Decimalised latitude for first point.
#' @param long2 Decimalised longitude for second point.
#' @param lat2 Decimalised latitude for second point.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return Great circle distance in kilometres.
#' @details Nothing yet.
#' @examples
#' # Return one-degree of latitude (approximately 111 kilometres):
#' GreatCircleDistanceFromLongLat(0, 0, 0, 1)

# GreatCircleDistance function:
GreatCircleDistanceFromLongLat <- function(long1, lat1, long2, lat2, EarthRad = 6367.4447) {
	
	# Convert point 1 longitude to radians:
	long1 <- long1 * (pi / 180)
	
	# Convert point 1 latitude to radians:
	lat1 <- lat1 * (pi / 180)
	
	# Convert point 2 longitude to radians:
	long2 <- long2 * (pi / 180)
	
	# Convert point 2 latitude to radians:
	lat2 <- lat2 * (pi / 180)
	
	# Get great circle distacne using cosine:
	great_circle_distance <- acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(abs(long1 - long2))) * EarthRad
	
	# Return great circle distance:
	return(great_circle_distance)
	
}
