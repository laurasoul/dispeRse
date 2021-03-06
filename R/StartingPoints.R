#' Generate a matrix of longitudes and latitudes of continental centres
#'
#' Generates a matrix of longitudes and latitudes describing the centres of circular continents on a sphere
#'
#' @param N_continents The (maximum) number of individual continents
#' @param radius The radius of each circular continent.
#' @param start_configuration One of "random separate", "random overlap", "supercontinent", or "max separate".
#' @param squishiness A value from 0 (continents can never overlap) to 1 (continents can overlap completely)
#' @param EarthRad Eartn radius in kilometres.
#' @param polar TRUE/FALSE Is there a continent starting on the south pole?
#' @return A matrix of longitudes and latitudes describing the centres of circular continents.
#' @details Nothing yet.
#' @author Laura C. Soul \email{lauracsoul@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' StartingPoints(N_continents = 7, radius = 2000,
#'    start_configuration = "supercontinent", squishiness = 0.1,
#'    EarthRad = 6367.4447, polar=FALSE)

StartingPoints <- function(N_continents = 7, radius = 2000, start_configuration = "supercontinent", squishiness = 0.25, EarthRad = 6367.4447, polar=FALSE) {

	# Check N_continents at least 1 and is an integer (divisble by one:
	if(N_continents %% 1 != 0 || N_continents < 1) stop("ERROR: Number of contientns must be a positive integer >= 1.")
	
	# Check radius is less than or equal to pi * the radius of the Earth:
	if(radius > (pi * EarthRad)) stop("ERROR: Radius must be less than pi * EarthRad (and much more if using multiple continents).")
	
	# Check configuration option chosen is valid:
	if(start_configuration != "random separate" && start_configuration != "random overlap" && start_configuration != "supercontinent" && start_configuration != "max separate") stop("ERROR: Starting configuration must be one of \"random separate\", \"random overlap\", \"supercontinent\", or \"max separate\".")
	
	# Check squishiness is a proportion:
	if(squishiness > 1 || squishiness < 0) stop("ERROR: Squishiness is proportional and must be between 0 and 1.")

	# If starting dispersed and the number of continents is larger that 8 use random distribution instead
	if(start_configuration == "max separate" && N_continents >= 9) {
	
		cat("Using random separate configuration as number of continents is too large", "\n")
		start_configuration <- "random separate"
	}
	
	# If the start configuration is a supercontinnet (all continents attached togther):
	if(start_configuration == "supercontinent") {
	
		# As long as squishiness is greater than 0 (complete overlap is impossible):
		if(squishiness > 0) {
		
			# Establish minimum distance between continent centres (:
			min_separation <- (1 - squishiness) * radius * 2
			
			# Randomly assign starting longitude:
			first_circle_long <- stats::runif(1, min = -180, max = 180)
			
			# Randomly assign starting latitude:
			first_circle_lat <- stats::runif(1, min = -90, max = 90)
			
			# Create matrix to store circles:
			circles <- matrix(c(first_circle_long, first_circle_lat), ncol=2, dimnames=list(c(), c("Longitude", "Latitude")))
			
			# If there are two or more continents still to place:
			if(N_continents > 1) {
				
				# Find centre of second continnet by randomly drawing from the circle describing all points within minimum separation of the first continent:
				second_circle <- EndPoint(longitude = first_circle_long, latitude = first_circle_lat, bearing = stats::runif(1, min = 0, max = 360), distance = min_separation, EarthRad = EarthRad)
				
				# Add second continent to circles matrix:
				circles <- rbind(circles, c(second_circle$longitude, second_circle$latitude))
				
				# If there are more than two continents still to place:
				if(N_continents > 2) {
					
					# Keep adding continents unless they all have starting points:
					while(N_continents > nrow(circles)) {
						
						# Find open spots available for next continent:
						open_spots <- OpenSpots(circles[, 1], circles[, 2], min_separation, EarthRad)
						
						# Error message if too many continents are requested:
						if(nrow(open_spots) == 0) stop("ERROR: There is not enough space for all your continents! Try making them smaller or choosing fewer of them.")

						# Randomly chose an open spot for the next continent and add to circles:
						circles <- rbind(circles, open_spots[sample(1:nrow(open_spots))[1], ])
						
					}
					
				}
				
			}
			
		# If squishiness is 1:
		} else {
			
			# Warn user:
			print("WARNING: Squishiness set to one and supercontinent option employed: all continents will start completely overlapping each other.")
			
			# Randomly assign starting longitude:
			first_circle_long <- stats::runif(1, min = -180, max = 180)
			
			# Randomly assign starting latitude:
			first_circle_lat <- stats::runif(1, min = -90, max = 90)
			
			# Create matrix to store circles:
			circles <- cbind(matrix(rep(c(first_circle_long, first_circle_lat), N_continents), ncol=2, byrow=TRUE))
			
			# Add column names:
			colnames(circles) <- c("Longitude", "Latitude")
			
		}
		
	}
	
	# If starting continental configuration is random:
	if(start_configuration == "random separate" || start_configuration == "random overlap") {
		
		# Establish minimum difference between continent centres if overlap allowed:
		if(start_configuration == "random overlap") min_separation <- (1 - squishiness) * radius * 2
		
		# Establish minimum difference between continent centres if overlap not allowed:
		if(start_configuration == "random separate") min_separation <- radius * 2
		
		# Error check for whether all continents can theoretically fit:
		if((SphericalCapArea(min_separation / 2) * N_continents) > (4 * pi * (EarthRad ^ 2))) stop("ERROR: The current choices for number of continents, radius and squishiness means it is impossible to fit all continents on the sphere. Consider choosing fewer continents, a smaller radius, or larger squishiness.")
		
		# Randomly assign starting longitude:
		first_circle_long <- stats::runif(1, min = -180, max = 180)
		
		# Randomly assign starting latitude:
		first_circle_lat <- stats::runif(1, min = -90, max = 90)
		
		# Create matrix to store circles:
		circles <- matrix(c(first_circle_long, first_circle_lat), ncol=2, dimnames=list(c(), c("Longitude", "Latitude")))
		
		# Keep adding continents until they all have starting points:
		while(N_continents > nrow(circles)) {

			# Randomly assign starting longitude:
			new_circle_long <- stats::runif(1, min = -180, max = 180)
			
			# Randomly assign starting latitude:
			new_circle_lat <- stats::runif(1, min = -90, max = 90)
			
			# Get the Great Circle distance matrix of the potential new continent and the pre-existing continents:
			GC_distance_matrix <- GreatCircleDistanceMatrix(rbind(circles[, c("Longitude", "Latitude")], c(new_circle_long, new_circle_lat))[, 1], rbind(circles[, c("Longitude", "Latitude")], c(new_circle_long, new_circle_lat))[, 2], EarthRad = EarthRad)
			
			# Reset counter:
			counter <- 1
			
			# As long as the new point is too close to the epre-existing continents:
			while(min(GC_distance_matrix[lower.tri(GC_distance_matrix)]) < min_separation) {
				
				# Randomly assign new starting longitude:
				new_circle_long <- stats::runif(1, min = -180, max = 180)
				
				# Randomly assign new starting latitude:
				new_circle_lat <- stats::runif(1, min = -90, max = 90)
				
				# Update Great Circle distance matrix:
				GC_distance_matrix <- GreatCircleDistanceMatrix(rbind(circles[, c("Longitude", "Latitude")], c(new_circle_long, new_circle_lat))[, 1], rbind(circles[, c("Longitude", "Latitude")], c(new_circle_long, new_circle_lat))[, 2], EarthRad = EarthRad)
				
				# Update counter:
				counter <- counter + 1
				
				# Break loop if it looks like adding the next continent is impossible:
				if(counter == 1000) stop("ERROR: Have made 1000 attempts to add the next continent without success. Check that number of continents or continent radius are not too large.")
				
			}
			
			# Add new continent to list:
			circles <- rbind(circles, c(new_circle_long, new_circle_lat))
			
		}
		
	}

	# If starting continental configuration is maximally separated:
	if(start_configuration == "max separate") {
		N <- N_continents
		#Puts the first continent on the south pole
		if (polar) {
			slong = 0
			slat = -90
		} else {
			#Randomly positions the first continent
			slong = stats::runif(1, c(-180,180))
			slat = stats::runif(1, c(-90, 90))
		}
		#matrix of continent coordinates
		points <- matrix(nrow=N, ncol=2)
		points[1,1] <- slong
		points[1,2] <- slat

		#Continents on either end of a diameter
		if (N==2)  {
			points[2,1] <- (slong - 180) %% 180
    		points[2,2] <- -slat
		}

		#Three points on a great circle
		if (N==3) {
			d<- 2*pi*EarthRad/3
			bearing <- stats::runif(1, c(0,360))
			second <- EndPoint(slong, slat, bearing, distance = d, EarthRad = EarthRad)
			third <- EndPoint(slong, slat, bearing, distance = 2 * d, EarthRad = EarthRad)
			points[2,] <- c(second$longitude, second$latitude)
			points[3,] <- c(third$longitude, third$latitude)
		}

		#Verticies of a tetrahedron
		if (N==4) {
			a<-EarthRad/(sqrt(3/8))
			theta <- ThetaFromChordLength(a)
			d <- EarthRad*theta
			bear <- stats::runif(1, c(0,360))
			second <- EndPoint(slong, slat, bearing=bear, distance = d, EarthRad = EarthRad)
			points[2,] <- c(second$longitude, second$latitude)
			third <- EndPoint(slong, slat, bearing=(bear+120) %% 360, distance = d, EarthRad = EarthRad)
			points[3,] <- c(third$longitude, third$latitude)
			fourth <- EndPoint(slong, slat, bearing=(bear+240) %% 360, distance = d, EarthRad = EarthRad)
			points[4,] <- c(fourth$longitude, fourth$latitude)
		}

		#No unique solution so three points on great circle and two at ends of perpendicular diameter
		if (N==5) {
			d<- 2*pi*EarthRad/3
			bear <- stats::runif(1, c(0,360))
			second <- EndPoint(slong, slat, bear, distance = d, EarthRad = EarthRad)
			third <- EndPoint(slong, slat, bear, distance = 2*d, EarthRad = EarthRad)
			fourth <- EndPoint(slong, slat, (bear + 90) %% 360, distance=0.5*EarthRad*pi, EarthRad = EarthRad)
			fifth <- c((fourth$longitude - 180) %% 180, -fourth$latitude)
			points[2,] <- c(second$longitude, second$latitude)
			points[3,] <- c(third$longitude, third$latitude)
   			points[4,] <- c(fourth$longitude, fourth$latitude)
    		points[5,] <- fifth
		}

		#Verticies of an octohedron
		if (N==6) {
			d<- 0.5 * pi * EarthRad
			bear <- stats::runif(1, c(0,360))
			second <- EndPoint(slong, slat, bear, distance = d, EarthRad = EarthRad)
			third <- EndPoint(slong, slat, bear, distance = 2*d, EarthRad = EarthRad)
			fourth <- EndPoint(slong, slat, bear, distance = 3*d, EarthRad = EarthRad)
			fifth <- EndPoint(slong, slat, (bear + 90) %% 360, distance=0.5*EarthRad*pi, EarthRad = EarthRad)
			sixth <- c((fifth$longitude - 180) %% 180, -fifth$latitude)
			points[2,] <- c(second$longitude, second$latitude)
			points[3,] <- c(third$longitude, third$latitude)
   			points[4,] <- c(fourth$longitude, fourth$latitude)
    		points[5,] <- c(fifth$longitude, fifth$latitude)
    		points[6,] <- sixth
		}
		#Verticies of a pentagonal bipyramid
		if (N==7) {
			d<- 2*pi*EarthRad/5
			bear <- stats::runif(1, c(0,360))
			second <- EndPoint(slong, slat, bear, distance = d, EarthRad = EarthRad)
			third <- EndPoint(slong, slat, bear, distance = 2*d, EarthRad = EarthRad)
			fourth <- EndPoint(slong, slat, bear, distance = 3*d, EarthRad = EarthRad)
			fifth <- EndPoint(slong, slat, bear, distance = 4*d, EarthRad = EarthRad)
			sixth <- EndPoint(slong, slat, (bear + 90) %% 360, distance=0.5*EarthRad*pi, EarthRad = EarthRad)
			seventh <- c((sixth$longitude - 180) %% 180, -sixth$latitude)
			points[2,] <- c(second$longitude, second$latitude)
			points[3,] <- c(third$longitude, third$latitude)
   			points[4,] <- c(fourth$longitude, fourth$latitude)
    		points[5,] <- c(fifth$longitude, fifth$latitude)
    		points[6,] <- c(sixth$longitude, sixth$latitude)
    		points[7,] <- seventh
		}
		# Verticies of a cube
		if (N==8) {
			d <- 2 * EarthRad * asin(1 / sqrt(3))
			second <- EndPoint(slong, slat, bearing=0, distance = d, EarthRad = EarthRad)
			third <- EndPoint(slong, slat, bearing=120, distance = d, EarthRad = EarthRad)
			fourth <- EndPoint(slong, slat, bearing=240, distance = d, EarthRad = EarthRad)
			reverse <- BearingBetweenTwoLongLatPoints(second$longitude, second$latitude, slong, slat)
			fifth <- EndPoint(second$longitude, second$latitude, bearing = (reverse + 120) %% 360, distance = d, EarthRad = EarthRad)
			sixth <- EndPoint(second$longitude, second$latitude, bearing = (reverse + 240) %% 360, distance = d, EarthRad = EarthRad)
			reverse2 <- BearingBetweenTwoLongLatPoints(sixth$longitude, sixth$latitude, second$longitude, second$latitude)
			seventh <- EndPoint(sixth$longitude, sixth$latitude, bearing = (reverse2 + 120) %% 360, distance = d, EarthRad = EarthRad)
			reverse3 <- BearingBetweenTwoLongLatPoints(third$longitude, third$latitude, slong, slat)
			eigth <- EndPoint(third$longitude, third$latitude, bearing = (reverse3 - 120) %% 360, distance = d, EarthRad = EarthRad)
			points[2,] <- c(second$longitude, second$latitude)
			points[3,] <- c(third$longitude, third$latitude)
   			points[4,] <- c(fourth$longitude, fourth$latitude)
    		points[5,] <- c(fifth$longitude, fifth$latitude)
    		points[6,] <- c(sixth$longitude, sixth$latitude)
    		points[7,] <- c(seventh$longitude, seventh$latitude)
    		points[8,] <- c(eigth$longitude, eigth$latitude)
		}

		# Create matrix to store circles:
		circles <- cbind(points)
			
		# Add column names:
		colnames(circles) <- c("Longitude", "Latitude")
	}
	
	# Output the circles matrix:
	return(circles)
	
}
