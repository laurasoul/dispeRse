#' Find open spots for next continent in a supercontinent
#'
#' Takes position of existing continents and finds the closest possible position(s) for a new continent
#' @param longitudes Decimalised longitudes of existing continents.
#' @param latitudes Decimalised longitudes of existing continents.
#' @param min_separation The minimum separation between continents in kilomeres.
#' @param EarthRad Earth radius in kilometres.
#' @return A two-column (longitude and latitude) matrix of open spots.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' OpenSpots(c(0, 1), c(0, 0), 111.132874953663, EarthRad = 6367.4447)

OpenSpots <- function(longitudes, latitudes, min_separation, EarthRad = 6367.4447) {
	
	# Get spherical angle of equilateral (spherical) triangle with sides equal to minimum separation:
	spherical_angle <- SphericalAngleForEquilateralTriangleFromGreatCircleSideLength(min_separation, EarthRad = EarthRad)
	
	# Get links between continents (essentially which ones are adjacent to each other):
	intercontinental_links <- IntercontinentalLinks(min_separation, longitudes, latitudes, EarthRad = EarthRad)
	
	# Create empty matrix to stoe open spots:
	open_spots <- matrix(nrow=0, ncol=2)
	
	# For each row in the lower triangle of intercontinental links:
	for(i in 2:nrow(intercontinental_links)) {
		
		# For each column in the lower triangle of intercontinental links:
		for(j in 1:(i - 1)) {
			
			# If there is a link:
			if(intercontinental_links[i, j] == 1) {
				
				# Get first of two possible bearings from centre of first continent to centre of adjacent continent:
				new_bearing_1 <- (BearingBetweenTwoLongLatPoints(longitudes[i], latitudes[i], longitudes[j], latitudes[j]) - spherical_angle) %% 360
				
				# Get second of two possible bearings from centre of first continent to centre of adjacent continent:
				new_bearing_2 <- (BearingBetweenTwoLongLatPoints(longitudes[i], latitudes[i], longitudes[j], latitudes[j]) + spherical_angle) %% 360
				
				# Get lat-long coordinates of first possible site new continent:
				new_continent_1 <- EndPoint(longitudes[i], latitudes[i], new_bearing_1, min_separation)[c("long", "lat")]
				
				# Get lat-long coordinates of second possible site new continent:
				new_continent_2 <- EndPoint(longitudes[i], latitudes[i], new_bearing_2, min_separation)[c("long", "lat")]
				
				# Add two new spots to open spots list:
				open_spots <- rbind(open_spots, rbind(unlist(new_continent_1), unlist(new_continent_2)))
				
			}
			
		}
		
	}
	
	# Ensure only unique open spots are included:
	open_spots <- matrix(as.numeric(unlist(strsplit(unique(apply(open_spots, 1, paste, collapse="%%")), "%%"))), ncol=2, byrow=TRUE)
	
	# Vector to store spots that are too close to (or already are occupied by) existing continents):
	too_close <- vector(mode="numeric")
	
	# For each open spot:
	for(i in 1:nrow(open_spots)) {
		
		# Get the distance from the ith open spot to all existing continents:
		distance_to_existing <- One2ManyGreatCircleDistance(open_spots[i, 1], open_spots[i, 2], longitudes, latitudes)
		
		# Isolate potential too short distances (need to check if just a floating point error before confirming):
		distance_to_existing <- distance_to_existing[distance_to_existing < min_separation]
		
		# If the ith open spot is potentially too close to an existing continent:
		if(length(distance_to_existing) > 0) {
			
			# For each potnetial too short distance::
			for(j in length(distance_to_existing):1) {
				
				# Exclude if within floating point error of the minimum allowed separation:
				if(all.equal(distance_to_existing[j], min_separation) == TRUE) distance_to_existing <- distance_to_existing[-j]
				
			}
			
		}
		
		# If the open spot is still too close to an existing continent add it to the list:
		if(length(distance_to_existing) > 0) too_close <- c(too_close, i)
		
	}
	
	# If any open spots are too close to existing continents then exclude them:
	if(length(too_close) > 0) open_spots <- open_spots[-too_close, ]
	
	# Return open spots matrix (may have zero rows):
	return(open_spots)
	
}
