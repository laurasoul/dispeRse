#' Which continents are linked?
#'
#' Returns a continent-by-continent matrix with intercontinental links scored as ones
#'
#' @param min_separation Minimum separation allowed between continental centres in kilometres.
#' @param longitudes A vector of the decimalised longitudes for each continental centre.
#' @param latitudes A vector of the decimalised latitudes for each continental centre.
#' @param EarthRad Earth radius in kilometres.
#' @return A square matrix with intercontinental links scored as ones.
#' @details Nothing yet.
#'
#' @export
#' @examples
#' IntercontinentalLinks(0, c(179, -179), c(-89, 89))

IntercontinentalLinks <- function(min_separation, longitudes, latitudes, EarthRad = 6367.4447) {
	
	# Get starting intercontinental distances:
	intercontinental_links <- intercontinental_distance_matrix <- GreatCircleDistanceMatrix(longitudes, latitudes, EarthRad = EarthRad)
	
	# For each row:
	for(i in 1:(nrow(intercontinental_distance_matrix) - 1)) {
		
		# For each column:
		for(j in (i + 1):nrow(intercontinental_distance_matrix)) {
			
			# If continents are linked:
			if(all.equal(intercontinental_distance_matrix[i, j], min_separation) == TRUE) {
				
				# Record a one for the link:
				intercontinental_links[j, i] <- intercontinental_links[i, j] <- 1
				
			# If continents are not linked:
			} else {
				
				# Record a zero for no link:
				intercontinental_links[j, i] <- intercontinental_links[i, j] <- 0
				
			}
			
		}
		
	}
	
	# Return intercontinental links matrix:
	return(intercontinental_links)
	
}
