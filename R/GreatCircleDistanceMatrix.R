#' Great Circle Distance Matrix Between Coordinates
#'
#' Returns a distance matrix between points of longitude and latitude
#' @param longs A vector of decimalised longitudes.
#' @param lats A vector of decimalised latitudes.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return A distancematrix of minimum Great Circle Distances.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' GreatCircleDistanceMatrix(c(-179, 0, 179), c(89, 0, -89))

GreatCircleDistanceMatrix <- function(longs, lats, EarthRad = 6367.4447) {
	
	# Create empty square matrix:
	GC_dist_matrix <- matrix(0,nrow = length(longs), ncol = length(longs))
	
	# For each first coordinate in pariwise comparison:
	for(i in 1:(length(longs) - 1)) {
		
		# For each second coordinate in pariwise comparison:
		for(j in (i + 1):length(longs)) {
			
			# Store minimal Great Circle distance:
			GC_dist_matrix[i, j] <- GC_dist_matrix[j, i] <- GreatCircleDistanceFromLongLat(longs[i], lats[i], longs[j], lats[j], EarthRad = EarthRad, Warn = FALSE)
			
		}
		
	}
	
	# Return Great Circle Distacne matrix with diagonal left as zero:
	return(GC_dist_matrix)
	
}
