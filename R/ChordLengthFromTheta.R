#' Returns chord length from theta
#'
#' Returns chord length in kilometres from theta in radians
#' @param theta Theta in radians.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return Chord length in kilometres.
#' @details Nothing yet.
#' @export
#' @examples
#' ChordLengthFromTheta(0.5)

# Get chord length from theta in radians and radius in km:
ChordLengthFromTheta <- function(theta, EarthRad = 6367.4447) {
	
	# Get chord length from theta (in radians):
	sin(theta / 2) * EarthRad * 2
	
}
