#' Equal-area grid from lines of longitude and latitude
#'
#' Makes an equal-area grid of specified size from lines of longitude and latitude
#'
#' @param N_longitude Theta in radians.
#' @param N_latitude Radius of the Earth in kilometres.
#' @return The lines of longitudde and latitude that describe the grid.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' EqualAreaRectangularGrid(N_latitude = 10, N_longitude = 20)
#' # map results
#' # Mercator projection

EqualAreaRectangularGrid <- function(N_latitude, N_longitude) {
	
	# Longitude boundaries:
	lonbreaks <- c(180, 180 - cumsum(rep(360 / N_longitude, N_longitude)))
	
	# Latitude boundaries:
	latbreaks <- asin(seq(-1, 1, length.out = N_latitude + 1)) * (180 / pi)
	
	# Combine results:
	result <- list(latbreaks, lonbreaks)
	
	# Add names:
	names(result) <- c("latitude_breaks", "longitude_breaks")
	
	# Return output:
	return(result)
	
}
