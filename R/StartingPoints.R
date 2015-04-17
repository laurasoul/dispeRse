#' Generate a matrix of longitudes and latitudes of contienntal centres
#'
#' Generates a matrix of longitudes and latitudes describing the centres of circular continents on a sphere
#'
#' @param N_continents The (maximum) number of individual continents
#' @param radius The radius of each circular continent.
#' @param start_configuration One of "random", "supercontinent", or "max separate".
#' @param squishiness A value from 0 (continents can never overlap) to 1 (continents can overlap completely)
#' @return A matrix of longitudes and latitudes describing the centres of circular continents
#' @details Nothing yet.
#' @examples
#' StartingPoints(N_continents = 3, radius = 1000, start_configuration = "supercontinent", squishiness = 0.1)

# Inouts for eventual continental function:
# - N circles
# - Radius (single value)
# - Starting points (random; supercontinent; all separate)
# - Squishiness (How close can they get?; %age of radii; make sure it is 0-100)
# - Starting bearing (random; disperse; converging) - will need to use Euler Pole function to do this.
# - Speed - Not sure how we will do this.
# - Stickiness: how long (in steps) can continents stay at maximum squished-togetherness?

# Outputs:
# - N separated continents (distance matrix with values less than 2radii)
# - Long-lat of each circle centre
# - Bearings after each step change
# - Total land area all circles - minus

StartingPoints <- function(N_continents = 3, radius = 500, start_configuration = "supercontinent", squishiness = 0.25) {

# Check N_continents at least 1 and is integer
# Check start_configuration is one of: "random", "supercontinent" "all separate"
# Check radius is positive and less than area of planet! (check with N_continents for total size)
# Check squishiness is between 0 and 1 (0 being 0%, circles cannot occlude at all; 1 being 100% circles can completely occlude)

	# If the start configuration is a supercontinnet (all continents attached togther):
	if(start_configuration == "supercontinent") {
	
		# Establish minimum distance between continent centres (:
		min_separation <- (1 - squishiness) * radius * 2
		
# If squishiness is zero then all continents begin with same centre!!!!!!!
# Need conditional to reflect this.
		
		# Randomly assign starting longitude:
		first_circle_long <- runif(1, min = -180, max = 180)

		# Randomly assign starting latitude:
		first_circle_lat <- runif(1, min = -90, max = 90)
		
		# Create matrix to store circles:
		circles <- matrix(c(1, first_circle_long, first_circle_lat), ncol=3, dimnames=list(c(), c("Circle", "Longitude", "Latitude")))
		
		# If there are two or more continents still to place:
		if(N_continents > 1) {
			
			# Find centre of second continnet by randomly drawing from the circle describing all points within minimum separation of the first continent:
			second_circle <- EndPoint(slong = first_circle_long, slat = first_circle_lat, bearing = runif(1, min = 0, max = 360), distance = min_separation)

			# Add second continent to circles matrix:
			circles <- rbind(circles, c(2, second_circle$long, second_circle$lat))
			
			# If there are more than two continents still to place:
			if(N_continents > 2) {
				
				# Get spherical angle of equilateral (spherical) triangle with sides equal to minimum separation:
				spherical_angle <- GetSphericalAngleForEquilateralTriangleFromGreatCircleSideLength(min_separation)

				# Get first of two possible bearings from centre of first continent to centre of third continent:
				new_bearing_1 <- (GetBearingBetweenTwoLongLatPoints(circles[1, "Longitude"], circles[1, "Latitude"], circles[2, "Longitude"], circles[2, "Latitude"]) - spherical_angle) %% 360
				
				# Get second of two possible bearings from centre of first continent to centre of third continent:
				new_bearing_2 <- (GetBearingBetweenTwoLongLatPoints(circles[1, "Longitude"], circles[1, "Latitude"], circles[2, "Longitude"], circles[2, "Latitude"]) + spherical_angle) %% 360
				
				# Get lat-long coordinates of first possible site of third continent:
				new_continent_1 <- EndPoint(circles[1, "Longitude"], circles[1, "Latitude"], new_bearing_1, min_separation)[c("long", "lat")]
				
				# Get lat-long coordinates of second possible site of third continent:
				new_continent_2 <- EndPoint(circles[1, "Longitude"], circles[1, "Latitude"], new_bearing_2, min_separation)[c("long", "lat")]
				
				# Randomly pick site for third continent:
				third_circle <- rbind(unlist(new_continent_1), unlist(new_continent_2))[sample(c(1, 2), 1), ]
				
				# Add third continent to circles matrix:
				circles <- rbind(circles, c(3, third_circle))
				
# Little plot to check things look like they are working OK!:
#library(maps)
#map()
#points(circles[, "Longitude"], circles[, "Latitude"], pch=19, col="red")
				
				# If there are more than three continents:
				if(N_continents > 3) {
				
					# Keep adding continents unless they all have starting points:
					while(N_continents > nrow(circles)) {
						
						
						
					}

# Now is a triangle problem only
# Establish centre of triangle first?
					
				}
				
			# If only two continents:
			} else {
			
# Output
				
			}
		
		# If only one continent:
		} else {
			
# Output
			
		}
		
	}
	
	# If starting continental configuration is random:
	if(start_configuration == "random") {
		
		# Establish minimum difference between points:
		min_separation <- (1 - squishiness) * radius * 2
		
	}

	# If starting continental configuration is maximally separated:
	if(start_configuration == "max separate") {
		
		# Establish minimum difference between points:
		#min_separation <- ?????
		
	}
	
	# Output the circles matrix:
	return(circles)
	
}
