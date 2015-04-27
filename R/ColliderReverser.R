#' Reverses a time step to coliision point
#'
#' If a collision happens within a time bin establishes the proportion of time elapsed
#'
#' @param min_separation The minimum separation between continents in kilometres.
#' @param continent_1_longitude_t0 Decimalised longitude of first continent at t0.
#' @param continent_1_latitude_t0 Decimalised latitude of first continent at t0.
#' @param continent_1_longitude_t1 Decimalised longitude of first continent at t1.
#' @param continent_1_latitude_t1 Decimalised latitude of first continent at t1.
#' @param continent_2_longitude_t0 Decimalised longitude of second continent at t0.
#' @param continent_2_latitude_t0 Decimalised latitude of second continent at t0.
#' @param continent_2_longitude_t1 Decimalised longitude of second continent at t1.
#' @param continent_2_latitude_t1 Decimalised latitude of second continent at t1.
#' @param continent_1_euler_longitude Decimalised longitude of Euler pole for first continent.
#' @param continent_1_euler_latitude Decimalised latitude of Euler pole for first continent.
#' @param continent_2_euler_longitude Decimalised longitude of Euler pole for second continent.
#' @param continent_2_euler_latitude Decimalised latitude of Euler pole for second continent.
#' @param continent_1_degrees_per_step Degrees per step (speed) of first continent.
#' @param continent_2_degrees_per_step Degrees per step (speed) of second continent.
#' @param EarthRad Radius of the Earth in kilometres.
#' @param Warn Whether or not to print warnings.
#' @return Proportion (0 to 1) of time step at which the minimum separation distance collision occurs.
#' @details Nothing yet.
#' @examples
#' min_separation <- 500
#'
#' continent_1_euler_longitude <- 2
#' continent_1_euler_latitude <- 89
#'
#' continent_2_euler_longitude <- 4
#' continent_2_euler_latitude <- -89
#'
#' continent_1_degrees_per_step <- 1
#' continent_2_degrees_per_step <- 1
#'
#' continent_1_longitude_t0 <- 3
#' continent_1_latitude_t0 <- 0
#' continent_2_longitude_t0 <- -2
#' continent_2_latitude_t0 <- 1
#'
#' continent_1_start_bearing <- 
#'   BearingBetweenTwoLongLatPoints(continent_1_euler_longitude,
#'   continent_1_euler_latitude, continent_1_longitude_t0, continent_1_latitude_t0)
#' continent_2_start_bearing <-
#'   BearingBetweenTwoLongLatPoints(continent_2_euler_longitude,
#'   continent_2_euler_latitude, continent_2_longitude_t0, continent_2_latitude_t0)
#' continent_1_euler_distance <-
#'   GreatCircleDistanceFromLongLat(continent_1_euler_longitude, continent_1_euler_latitude,
#'   continent_1_longitude_t0, continent_1_latitude_t0)
#' continent_2_euler_distance <-
#'   GreatCircleDistanceFromLongLat(continent_2_euler_longitude, continent_2_euler_latitude,
#'   continent_2_longitude_t0, continent_2_latitude_t0)
#'
#' continent_1_longitude_t1 <- EndPoint(continent_1_euler_longitude, continent_1_euler_latitude,
#'   continent_1_start_bearing + continent_1_degrees_per_step, continent_1_euler_distance)$long
#' continent_1_latitude_t1 <- EndPoint(continent_1_euler_longitude, continent_1_euler_latitude,
#'   continent_1_start_bearing + continent_1_degrees_per_step, continent_1_euler_distance)$lat
#' continent_2_longitude_t1 <- EndPoint(continent_2_euler_longitude, continent_2_euler_latitude,
#'   continent_2_start_bearing + continent_2_degrees_per_step, continent_2_euler_distance)$long
#' continent_2_latitude_t1 <- EndPoint(continent_2_euler_longitude, continent_2_euler_latitude,
#'   continent_2_start_bearing + continent_2_degrees_per_step, continent_2_euler_distance)$lat
#'
#' GreatCircleDistanceFromLongLat(continent_1_longitude_t0, continent_1_latitude_t0,
#'   continent_2_longitude_t0, continent_2_latitude_t0)
#' GreatCircleDistanceFromLongLat(continent_1_longitude_t1, continent_1_latitude_t1,
#'   continent_2_longitude_t1, continent_2_latitude_t1)
#'
#' ColliderReverser(min_separation, continent_1_longitude_t0, continent_1_latitude_t0,
#'   continent_1_longitude_t1, continent_1_latitude_t1, continent_2_longitude_t0,
#'   continent_2_latitude_t0, continent_2_longitude_t1, continent_2_latitude_t1,
#'   continent_1_euler_longitude, continent_1_euler_latitude, continent_2_euler_longitude,
#'   continent_2_euler_latitude, continent_1_degrees_per_step, continent_2_degrees_per_step)

continent_1_longitude_t0 <- position[cont_1, t-1, 1]
continent_1_latitude_t0 <- position[cont_1, t-1, 2]
continent_1_longitude_t1 <- temp_position[cont_1,1]
continent_1_latitude_t1 <- temp_position[cont_1,2]
continent_2_longitude_t0 <- position[cont_2, t-1, 1]
continent_2_latitude_t0 <- position[cont_2, t-1, 2]
continent_2_longitude_t1 <- temp_position[cont_2,1]
continent_2_latitude_t1 <- temp_position[cont_2,2]


ColliderReverser <- function(min_separation, continent_1_longitude_t0, continent_1_latitude_t0, continent_1_longitude_t1, continent_1_latitude_t1, continent_2_longitude_t0, continent_2_latitude_t0, continent_2_longitude_t1, continent_2_latitude_t1, continent_1_euler_longitude, continent_1_euler_latitude, continent_2_euler_longitude, continent_2_euler_latitude, continent_1_degrees_per_step, continent_2_degrees_per_step, EarthRad = 6367.4447, Warn = TRUE) {
	
	# Get t0 distance between continents for checks on input data:
	t0_distance <- GreatCircleDistanceFromLongLat(continent_1_longitude_t0, continent_1_latitude_t0, continent_2_longitude_t0, continent_2_latitude_t0, Warn = FALSE)
	
	# Get t1 distance between continents for checks on input data:
	t1_distance <- GreatCircleDistanceFromLongLat(continent_1_longitude_t1, continent_1_latitude_t1, continent_2_longitude_t1, continent_2_latitude_t1, Warn = FALSE)

	# Check that continents did not converge before time window:
	if(min_separation > t0_distance) stop("Continents convrged before t0.")
	
	# Check that continents do converge in the interval:
	if(t0_distance < t1_distance) stop("Continents are not converging between t0 and t1.")

	# Check that continents cross minimum separation threshold in time window (i.e., between t0 and t1):
	if(min_separation < t1_distance) stop("Continents do not cross minimum separation threshold within time window")
	
	# If minimum separation threshold has been reached at start (t0):
	if(all.equal(t0_distance, min_separation)[1] == TRUE) {
		
		# Set degree modifier to zero:
		degree_modifier <- 0
		
		# Warn user of outcome:
		if(Warn) print("WARNING: Continents have already collided (within floating point error) at t0. Returning zero.")

	}
	
	# If minimum separation threshold is reached at end (t1):
	if(all.equal(t1_distance, min_separation)[1] == TRUE) {
		
		# Set degree modifier to zero:
		degree_modifier <- 1
		
		# Warn user of outcome:
		if(Warn) print("WARNING: Continents collided (within floating point error) at t1. Returning one.")
		
	}
	
	# If minimum separation threshold is reached *between* t0 and t1:
	if(!exists("degree_modifier")) {
		
		# Get bearing from Euler pole for first continent to t0 position of first continent:
		continent_1_start_bearing <- BearingBetweenTwoLongLatPoints(continent_1_euler_longitude, continent_1_euler_latitude, continent_1_longitude_t0, continent_1_latitude_t0)
		
		# Get bearing from Euler pole for second continent to t0 position of second continent:
		continent_2_start_bearing <- BearingBetweenTwoLongLatPoints(continent_2_euler_longitude, continent_2_euler_latitude, continent_2_longitude_t0, continent_2_latitude_t0)
		
		# Get GC distance from Euler pole of first continent to first continent centre:
		continent_1_euler_distance <- GreatCircleDistanceFromLongLat(continent_1_euler_longitude, continent_1_euler_latitude, continent_1_longitude_t0, continent_1_latitude_t0, EarthRad = EarthRad, Warn = FALSE)
		
		# Get GC distance from Euler pole of second continent to second continent centre:
		continent_2_euler_distance <- GreatCircleDistanceFromLongLat(continent_2_euler_longitude, continent_2_euler_latitude, continent_2_longitude_t0, continent_2_latitude_t0, EarthRad = EarthRad, Warn = FALSE)
		
		# Scalar which describes the proportion of a step at which continents are exactly minimum separation apart (will be modified!):
		degree_modifier <- 0.5
		
		# Starting stepsize (will shrink as answer is honed in on more precisely):
		stepsize <- 0.1
		
		# Starting coordinate of degree-modified guess for the point at which the two continents are exactly minimally separated for first continent:
		x <- EndPoint(continent_1_euler_longitude, continent_1_euler_latitude, (continent_1_start_bearing + (continent_1_degrees_per_step * degree_modifier)) %% 360, continent_1_euler_distance)[c("long", "lat")]
		
		# Starting coordinate of degree-modified guess for the point at which the two continents are exactly minimally separated for second continent:
		y <- EndPoint(continent_2_euler_longitude, continent_2_euler_latitude, (continent_2_start_bearing + (continent_2_degrees_per_step * degree_modifier)) %% 360, continent_2_euler_distance)[c("long", "lat")]
		
		# As long as we have not reached the point where the distacne between the two continents 
		while(!all.equal(GreatCircleDistanceFromLongLat(x$long, x$lat, y$long, y$lat, EarthRad = EarthRad, Warn = FALSE), min_separation) == TRUE) {
			
			# First establish current distance between continents based on best guess for degree modifier:
			current_distance <- GreatCircleDistanceFromLongLat(x$long, x$lat, y$long, y$lat, Warn = FALSE)
			
			# Do we need to increase the degree modifier value?:
			if(GreatCircleDistanceFromLongLat(x$long, x$lat, y$long, y$lat, EarthRad = EarthRad, Warn = FALSE) > min_separation) {
				
				# Create potential better value by adding step size to modifier:
				limit <- degree_modifier + stepsize
				
			# Or do we need to decrease the degree modifier value?:
			} else {
				
				# Create potential better value by subtracting step size to modifier:
				limit <- degree_modifier - stepsize
				
			}
			
			# Get new x based on potential better modifer:
			new_x <- EndPoint(continent_1_euler_longitude, continent_1_euler_latitude, (continent_1_start_bearing + (continent_1_degrees_per_step * limit)) %% 360, continent_1_euler_distance)[c("long", "lat")]
			
			# Get new y based on potential better modifer:
			new_y <- EndPoint(continent_2_euler_longitude, continent_2_euler_latitude, (continent_2_start_bearing + (continent_2_degrees_per_step * limit)) %% 360, continent_2_euler_distance)[c("long", "lat")]
			
			# Get new distance based on potential better modifer:
			new_distance <- GreatCircleDistanceFromLongLat(new_x$long, new_x$lat, new_y$long, new_y$lat, EarthRad = EarthRad, Warn = FALSE)
			
			# If new distance is closer to minimum separation:
			if(abs(current_distance - min_separation) > abs(new_distance - min_separation)) {
				
				# Update x based on new modifier:
				x <- new_x
				
				# Update y based on new modifier:
				y <- new_y
				
				# Update degree modifier itself:
				degree_modifier <- limit
				
			# If current distacne is stil our best estimate:
			} else {
				
				# Shrink the step size so we can hone in closer:
				stepsize <- stepsize * 0.1
				
			}
			
		}
		
	}
	
	# Output degree modifier (proportion of time step needed to travel):
	return(degree_modifier)
	
}

