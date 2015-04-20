#' Generate a matrix of longitudes and latitudes of contienntal centres
#'
#' Generates a matrix of longitudes and latitudes describing the centres of circular continents on a sphere
#'
#' @param N_continents The (maximum) number of individual continents
#' @param radius The radius of each circular continent.
#' @param start_configuration One of "random separate", "random overlap", "supercontinent", or "max separate".
#' @param squishiness A value from 0 (continents can never overlap) to 1 (continents can overlap completely)
#' @return A matrix of longitudes and latitudes describing the centres of circular continents
#' @details Nothing yet.
#' @examples
#' StartingPoints(N_continents = 7, radius = 2000,
#'    start_configuration = "supercontinent", squishiness = 0.1,
#'    EarthRad = 6367.4447)

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

StartingPoints <- function(N_continents = 7, radius = 2000, start_configuration = "supercontinent", squishiness = 0.25, EarthRad = 6367.4447) {

	# Check N_continents at least 1 and is an integer (divisble by one:
	if(N_continents %% 1 != 0 || N_continents < 1) stop("ERROR: Number of contientns must be a positive integer >= 1.")
	
	# Check radius is less than or equal to pi * the radius of the Earth:
	if(radius > (pi * EarthRad)) stop("ERROR: Radius must be less than pi * EarthRad (and much more if using multiple continents).")
	
	# Check configuration option chosen is valid:
	if(start_configuration != "random separate" && start_configuration != "random overlap" && start_configuration != "supercontinent" && start_configuration != "all separate") stop("ERROR: Starting configuration must be one of \"random separate\", \"random overlap\", \"supercontinent\", or \"all separate\".")
	
	# Check squishiness is a proportion:
	if(squishiness > 1 || squishiness < 0) stop("ERROR: Squishiness is proportional and must be between 0 and 1.")

	# If the start configuration is a supercontinnet (all continents attached togther):
	if(start_configuration == "supercontinent") {
	
		# As long as squishiness is greater than 0 (complete overlap is impossible):
		if(squishiness > 0) {
		
			# Establish minimum distance between continent centres (:
			min_separation <- (1 - squishiness) * radius * 2
			
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
					spherical_angle <- SphericalAngleForEquilateralTriangleFromGreatCircleSideLength(min_separation, EarthRad = EarthRad)
					
					# Get first of two possible bearings from centre of first continent to centre of third continent:
					new_bearing_1 <- (BearingBetweenTwoLongLatPoints(circles[1, "Longitude"], circles[1, "Latitude"], circles[2, "Longitude"], circles[2, "Latitude"]) - spherical_angle) %% 360
					
					# Get second of two possible bearings from centre of first continent to centre of third continent:
					new_bearing_2 <- (BearingBetweenTwoLongLatPoints(circles[1, "Longitude"], circles[1, "Latitude"], circles[2, "Longitude"], circles[2, "Latitude"]) + spherical_angle) %% 360
					
					# Get lat-long coordinates of first possible site of third continent:
					new_continent_1 <- EndPoint(circles[1, "Longitude"], circles[1, "Latitude"], new_bearing_1, min_separation)[c("long", "lat")]
					
					# Get lat-long coordinates of second possible site of third continent:
					new_continent_2 <- EndPoint(circles[1, "Longitude"], circles[1, "Latitude"], new_bearing_2, min_separation)[c("long", "lat")]
					
					# Randomly pick site for third continent:
					third_continent_picked <- sample(c(1, 2), 1)
					
					# Randomly pick site for third continent:
					third_circle <- rbind(unlist(new_continent_1), unlist(new_continent_2))[third_continent_picked, ]
					
					# Add third continent to circles matrix:
					circles <- rbind(circles, c(3, third_circle))
					
					# If there are more than three continents:
					if(N_continents > 3) {
						
						# Create new matrix for open spots (centres where new continents may reside) with unused spot from third continent choice:
						open_spots <- matrix(rbind(unlist(new_continent_1), unlist(new_continent_2))[setdiff(c(1:2), third_continent_picked), ], ncol=2)
						
						# Add column names:
						colnames(open_spots) <- c("Longitude", "Latitude")
						
						# Vector to store distances from centres:
						distances_from_centres <- vector(mode="numeric")
						
						# For each pre-existing continent:
						for(i in 1:3) {
							
							# Get distance from open spot to ith continent centre:
							distances_from_centres[i] <- GreatCircleDistanceFromLongLat(open_spots[, "Longitude"], open_spots[, "Latitude"], circles[i, "Longitude"], circles[i, "Latitude"], EarthRad = EarthRad)
							
						}
						
						# Identify continent centre furthest from current single open spot:
						furthest_point <- match(max(distances_from_centres), distances_from_centres)
						
						# Get bearing from furthest point to first of the two other continent centres:
						first_bearing <- BearingBetweenTwoLongLatPoints(circles[furthest_point, "Longitude"], circles[furthest_point, "Latitude"], circles[setdiff(c(1:3), furthest_point)[1], "Longitude"], circles[setdiff(c(1:3), furthest_point)[1], "Latitude"])
						
						# Get bearing from furthest point to second of the two other continent centres:
						second_bearing <- BearingBetweenTwoLongLatPoints(circles[furthest_point, "Longitude"], circles[furthest_point, "Latitude"], circles[setdiff(c(1:3), furthest_point)[2], "Longitude"], circles[setdiff(c(1:3), furthest_point)[2], "Latitude"])
						
						# Get bearings from furthest point to the two new open spots:
						new_bearings <- as.vector(sort(c(first_bearing - spherical_angle, first_bearing + spherical_angle, second_bearing - spherical_angle, second_bearing + spherical_angle))[c(1, 4)] %% 360)
						
						# Add first new open spot to list:
						open_spots <- rbind(open_spots, as.vector(unlist(EndPoint(circles[furthest_point, "Longitude"], circles[furthest_point, "Latitude"], new_bearings[1], min_separation)[c("long", "lat")])))
						
						# Add second new open spot to list:
						open_spots <- rbind(open_spots, as.vector(unlist(EndPoint(circles[furthest_point, "Longitude"], circles[furthest_point, "Latitude"], new_bearings[2], min_separation)[c("long", "lat")])))
						
						# Keep adding continents unless they all have starting points:
						while(N_continents > nrow(circles)) {
							
							# Error message if too many continents are requested:
							if(nrow(open_spots) == 0) stop("ERROR: There is not enough space for all your continents! Try making them smaller or choosing fewer of them.")
							
							# Randomly pick an open spot to place a new continent:
							new_continent_point <- sample(c(1:nrow(open_spots)))[1]
							
							# Vector to store row numbers for adjacent continents:
							adjacent_continents <- vector(mode="numeric")
							
							# For each pre-existing continent:
							for(i in 1:nrow(circles)) {
								
								# Is the distance from the new continent to the ith pre-existing continent equal to the minimum continental separation?:
								equal_check <- all.equal(as.vector(GreatCircleDistanceFromLongLat(open_spots[new_continent_point, "Longitude"], open_spots[new_continent_point, "Latitude"], circles[i, "Longitude"], circles[i, "Latitude"], EarthRad = EarthRad)), min_separation)
								
								# If one of the nearest continents then add to the list:
								if(equal_check == TRUE) adjacent_continents <- c(adjacent_continents, i)
								
							}
							
							# Vector to store bearings from new continent to adjacent continents:
							bearings_to_adjacents <- vector(mode="numeric")
							
							# For each continent adjacent to the new continent:
							for(i in adjacent_continents) {
								
								# Store bearing from new continent to ith adjacent continent:
								bearings_to_adjacents <- as.vector(c(bearings_to_adjacents, BearingBetweenTwoLongLatPoints(open_spots[new_continent_point, "Longitude"], open_spots[new_continent_point, "Latitude"], circles[i, "Longitude"], circles[i, "Latitude"])))
								
							}
							
							# Vector to store bearings to potential new open spots:
							new_bearings <- vector(mode="numeric")
							
							# For each bearing to an adjacent continent:
							for(i in 1:length(bearings_to_adjacents)) {
								
								# Get new bearings from plus or minus the psherical angle from the current ebarings:
								new_bearings <- c(new_bearings, c(bearings_to_adjacents[i] + spherical_angle, bearings_to_adjacents[i] - spherical_angle))
								
							}
							
							# Get bearings from new continent to two new potential open spots:
							new_bearings <- sort(new_bearings)[c(1, length(new_bearings))] %% 360
							
							# Store potential new open spots for further vetting:
							potential_new_open_spots <- rbind(as.vector(unlist(EndPoint(open_spots[new_continent_point, "Longitude"], open_spots[new_continent_point, "Latitude"], new_bearings[1], min_separation)[c("long", "lat")])), as.vector(unlist(EndPoint(open_spots[new_continent_point, "Longitude"], open_spots[new_continent_point, "Latitude"], new_bearings[2], min_separation)[c("long", "lat")])))
							
							# Add column names:
							colnames(potential_new_open_spots) <- c("Longitude", "Latitude")
							
							# Add new continent to circles matrix:
							circles <- rbind(circles, c(nrow(circles) + 1, open_spots[new_continent_point, ]))
							
							# Vector to store unsuitable new spots (those that are closer than minimum separation to existing continents):
							unsuitable_spots <- vector(mode="numeric")
							
							# For each potential new spot:
							for(i in 1:2) {
								
								# For each existing continent:
								for(j in 1:nrow(circles)) {
									
									# Get distance from ith potential new open spot to jth existing continent:
									distance_between_centres <- as.vector(GreatCircleDistanceFromLongLat(potential_new_open_spots[i, "Longitude"], potential_new_open_spots[i, "Latitude"], circles[j, "Longitude"], circles[j, "Latitude"], EarthRad = EarthRad))
									
									# If ith potential new open spot is closer to a pre-existing continent than the minimum separation allows add to unsuitable spots:
									if(all.equal(distance_between_centres, min_separation) != TRUE && distance_between_centres < min_separation) unsuitable_spots <- as.vector(c(unsuitable_spots, i))
									
								}
								
								# For each existing open spot:
								for(j in 1:nrow(open_spots)) {
									
									# Get distance from ith potential new open spot to jth existing open spot:
									distance_between_centres <- as.vector(GreatCircleDistanceFromLongLat(potential_new_open_spots[i, "Longitude"], potential_new_open_spots[i, "Latitude"], open_spots[j, "Longitude"], open_spots[j, "Latitude"], EarthRad = EarthRad))
									
									# If ith potential new open spot is closer to a pre-existing open spot than the minimum separation allows add to unsuitable spots:
									if(all.equal(distance_between_centres, min_separation) != TRUE && distance_between_centres < min_separation) unsuitable_spots <- as.vector(c(unsuitable_spots, i))
									
								}
								
							}
							
							# If there are unsuitable spots:
							if(length(unsuitable_spots) > 0) {
								
								# Remove unsuitable spots from potential new spots:
								potential_new_open_spots <- potential_new_open_spots[-unique(unsuitable_spots), ]
								
								# If deletions has de-matrixed the data:
								if(!is.matrix(potential_new_open_spots)) {
									
									# Make back into matrix:
									potential_new_open_spots <- matrix(potential_new_open_spots, ncol=2)
									
								}
								
								# Add column names:
								colnames(potential_new_open_spots) <- c("Longitude", "Latitude")
								
							}
							
							# Add new potential open spots to open spots:
							open_spots <- rbind(open_spots, potential_new_open_spots)
							
							# Remove new continent from open spots:
							open_spots <- open_spots[-new_continent_point, ]
							
							# If deletions has de-matrixed the data:
							if(!is.matrix(open_spots)) {
								
								# Make back into matrix:
								open_spots <- matrix(open_spots, ncol=2)
								
							}
							
							# Add column names:
							colnames(open_spots) <- c("Longitude", "Latitude")
							
						}
						
					}
					
				}
				
			}
			
		# If squishiness is 1:
		} else {
			
			# Warn user:
			print("WARNING: Squishiness set to one and supercontinent option employed: all continents will start completely overlapping each other.")
			
			# Randomly assign starting longitude:
			first_circle_long <- runif(1, min = -180, max = 180)
			
			# Randomly assign starting latitude:
			first_circle_lat <- runif(1, min = -90, max = 90)
			
			# Create matrix to store circles:
			circles <- cbind(c(1:N_continents), matrix(rep(c(first_circle_long, first_circle_lat), N_continents), ncol=2, byrow=TRUE))
			
			# Add column names:
			colnames(circles) <- c("Circle", "Longitude", "Latitude")
			
		}
		
	}
	
	# If starting continental configuration is random:
	if(start_configuration == "random separate" || start_configuration == "random overlap") {
		
		# Establish minimum difference between continent centres if overlap allowed:
		if(start_configuration == "random overlap") min_separation <- (1 - squishiness) * radius * 2
		
		# Establish minimum difference between continent centres if overlap not allowed:
		if(start_configuration == "random separate") min_separation <- radius * 2
		
		# Randomly assign starting longitude:
		first_circle_long <- runif(1, min = -180, max = 180)
		
		# Randomly assign starting latitude:
		first_circle_lat <- runif(1, min = -90, max = 90)
		
		# Create matrix to store circles:
		circles <- matrix(c(1, first_circle_long, first_circle_lat), ncol=3, dimnames=list(c(), c("Circle", "Longitude", "Latitude")))
		
		# Keep adding continents until they all have starting points:
		while(N_continents > nrow(circles)) {

			
			
			
		}
		
	}

	# If starting continental configuration is maximally separated:
	if(start_configuration == "max separate") {
		
		# Establish minimum difference between points:
		#min_separation <- ?????
		
	}
	
	# Output the circles matrix:
	return(circles)
	
}
