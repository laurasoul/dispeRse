#' Get theta from chord length
#'
#' Returns theta in radians from a chord length in kilometres
#'
#' @param chord_length The chord length in kilometres.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return Theta in radians.
#' @details Nothing yet.
#' @examples
#' ThetaFromChordLength(1000)

# Get theta in radians from chord length and radius in km:
ThetaFromChordLength <- function(chord_length, EarthRad = 6367.4447) {
	
	# Get theta in radians from chord length and radius in km:
	asin((chord_length * 0.5) / EarthRad) * 2
	
}
