#' Equal-area grid from lines of longitude and latitude
#'
#' Makes an equal-area grid of specified size from lines of longitude and latitude
#' @param N_longitude Number of East-West bins.
#' @param N_latitude Number of North-South bins.
#' @return The decimalised lines of longitude and latitude that describe the grid.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' # Create a 20-by-10 (200 cell) grid:
#' grids <- EqualAreaRectangularGrid(N_longitude = 20, N_latitude = 10)
#'
#' # View the output:
#' grids
#'
#' # Visualise the output using a Mercator projection:
#' plot(0, 0, type = "n", xlim = c(-180, 180), ylim = c(-90, 90), xlab = "Longitude", ylab = "Latitude")
#' for(i in 1:length(grids$latitude_breaks)) lines(c(-180, 180), c(grids$latitude_breaks[i], grids$latitude_breaks[i]), lty = 2)
#' for(i in 1:length(grids$longitude_breaks)) lines(c(grids$longitude_breaks[i], grids$longitude_breaks[i]), c(-90, 90), lty = 2)

EqualAreaRectangularGrid <- function(N_longitude, N_latitude) {
	
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
