#' Get theta from great circl distance
#'
#' Returns theta in radians from a great circle distance in kilometres
#'
#' @param great_circle_distance The great circle distance in kilometres.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return Theta in radians.
#' @details Nothing yet.
#' @examples
#' ThetaFromGreatCircleDistance(1000)

# Get theta in radians from great circle distance:
ThetaFromGreatCircleDistance <- function(great_circle_distance, EarthRad = 6367.4447) {
	
	# Get theta in radians from great circle distance:
	great_circle_distance / EarthRad
	
}
