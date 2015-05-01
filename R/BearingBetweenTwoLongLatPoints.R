#' Returns bearing from a set of coordinates
#'
#' Returns bearing in degrees for starting and ending longitude-latitude points
#' @param start_long Decimalised longitude for start point.
#' @param start_lat Decimalised latitude for start point.
#' @param end_long Decimalised longitude for end point.
#' @param end_lat Decimalised latitude for end point.
#' @return Bearing in degrees (0 or 360 for North, 90 for East, etc.).
#' @details Nothing yet.
#' @export
#' @examples
#' # Should return bearing of North (0):
#' BearingBetweenTwoLongLatPoints(start_long = 0, start_lat = 0,
#'   end_long = 0, end_lat = 1)

# Get bearing (with 0 or 360 being North, 90 being East etc.) from one lat-long coordinate to another:
BearingBetweenTwoLongLatPoints <- function(start_long, start_lat, end_long, end_lat) {
	
# Check that lat and long are not both identical (where no bearing is possible)
	
	# Convert starting longitude to radians:
	start_long <- start_long * (pi / 180)
	
	# Convert starting latitude to radians:
	start_lat <- start_lat * (pi / 180)
	
	# Convert ending longitude to radians:
	end_long <- end_long * (pi / 180)
	
	# Convert ending latitude to radians:
	end_lat <- end_lat * (pi / 180)
	
	# Get bearing in radians:
	bearing <- atan2(sin(end_long - start_long) * cos(end_lat), cos(start_lat) * sin(end_lat) - sin(start_lat) * cos(end_lat) * cos(end_long - start_long))
	
	# Convert bearing to degrees:
	bearing <- bearing / (pi / 180)
	
	# Place bearing on 0 to 360 scale, modulo 360:
	bearing <- (bearing + 360) %% 360
	
	# Return bearing:
	return(bearing)
	
}
