#' How many separate continents are there
#'
#' Given a minimum separation between continental centres records how many separate continents exist
#'
#' @param min_separation Minimum separation allowed between continental centres in kilometres.
#' @param longitudes A vector of the decimalised longitudes for each continental centre.
#' @param latitudes A vector of the decimalised latitudes for each continental centre.
#' @param EarthRad Earth radius in kilometres.
#' @return A vector of length equal to the number of separate continents.
#' @details Nothing yet.
#'
#' @export
#' @examples
#' HowManySeparateContinents(100, c(179, -179), c(-89, 89))

HowManySeparateContinents <- function(min_separation, longitudes, latitudes, EarthRad = 6367.4447) {
	
	# Get intercontinental links:
	intercontinental_links <- IntercontinentalLinks(min_separation, longitudes, latitudes, EarthRad = EarthRad) 
	
	# Add row and column names to matrix:
	colnames(intercontinental_links) <- rownames(intercontinental_links) <- 1:length(longitudes)
	
	# If there are no joined continents:
	if(sum(intercontinental_links) == 0) {
		
		# Simply output a vector of numbers (one for each separate continent):
		output <- colnames(intercontinental_links)
		
	# If there are joined continents:
	} else {
		
		# Create empty character vector to store each set of joined continents:
		output <- vector(mode="character")
		
		# Whilst there are remaining continental links to delineate:
		while(nrow(intercontinental_links) > 0) {
			
			# First find the (first) remaining continent with the most links to other continents:
			largest_clump <- colnames(intercontinental_links)[match(max(apply(intercontinental_links, 1, sum)), apply(intercontinental_links, 1, sum))]
			
			# If there are multiple continents define the clump by those continents that form direct links to the first continent:
			if(length(intercontinental_links) > 1) largest_clump <- sort(c(largest_clump, names(which(intercontinental_links[largest_clump, ] == 1))))
			
			# Enumerate remaining continents (not directly linked to first continent):
			remaining_continents <- setdiff(colnames(intercontinental_links), largest_clump)
			
			# If there are continents left to delineate:
			if(length(remaining_continents) > 0) {
				
				# For each remaining continent:
				for(i in remaining_continents) {
					
					# If the continent is linked to any continent already part of the current clump:
					if(length(sort(match(names(which(intercontinental_links[i, ] == 1)), largest_clump))) > 0) {
						
						# Add continent to clump:
						largest_clump <- c(largest_clump, i)
						
					}
					
				}
				
				# Update remaining continents:
				remaining_continents <- setdiff(colnames(intercontinental_links), largest_clump)
				
			}	
			
			# Add current continental clump to output:
			output <- c(output, paste(sort(largest_clump), collapse="&"))
			
			# Find rows to delete:
			delete_rows <- match(largest_clump, rownames(intercontinental_links))
			
			# Remove current continental clump from links matrix:
			intercontinental_links <- intercontinental_links[-delete_rows, -delete_rows]
			
			# Catch issue if no longer a matrix and force back into a 1-by-1 matrix:
			if(!is.matrix(intercontinental_links)) intercontinental_links <- matrix(intercontinental_links, ncol=1, dimnames=list(remaining_continents, remaining_continents))
			
		}
		
	}
	
	# Return unique continent list:
	return(output)
	
}
