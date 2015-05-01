#' Returns great circle distance from theta
#'
#' Returns great circle distance in kilometres from theta in radians
#' @param theta Theta in radians.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return Great circle distance in kilometres.
#' @details Nothing yet.
#' @export
#' @examples
#' GreatCircleDistanceFromTheta(0.5)

# Get great circle distacne from theta in radians and radius in km:
GreatCircleDistanceFromTheta <- function(theta, EarthRad = 6367.4447) {
	
	# Calculate great circle distance from theta in radians:
	theta * EarthRad
	
}
