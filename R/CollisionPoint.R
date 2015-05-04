#' Collision point
#'
#' Finds point at which an organism collides with the coast when dispersing
#'
#' @param start_longitude Decimalised longitude of organism's starting point.
#' @param start_latitude Decimalised latitude of organism's starting point.
#' @param end_longitude Decimalised longitude of organism's finishing point.
#' @param end_latitude Decimalised latitude of organism's finishing point.
#' @param continent_longitude Theta in radians.
#' @param continent_latitude Theta in radians.
#' @param continent_radius Continent radius in kilometres
#' @param EarthRad Radius of the Earth in kilometres.
#' @return Decimalised longitude and latitude of collision point and remaining unspent dispersal distance in kilometres.
#' @details Nothing yet. Circular continent.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' #CollisionPoint()

# Get chord length from theta in radians and radius in km:
CollisionPoint <- function(start_longitude, start_latitude, end_longitude, end_latitude, continent_longitude, continent_latitude, continent_radius, EarthRad = 6367.4447) {
	
# Check input data (does colllision actually occur?)
	
	# Scalar which describes the proportion of the distance the organism has to travel before colliding with the continent edge:
	distance_modifier <- 0.5
	
	# Starting stepsize (will shrink as answer is honed in on more precisely):
	stepsize <- 0.1
	
	# Get bearing of dispersal step:
	bearing <- BearingBetweenTwoLongLatPoints(start_longitude, start_latitude, end_longitude, end_latitude)

	# Get distance of dispersal step:
	distance <- GreatCircleDistanceFromLongLat(start_longitude, start_latitude, end_longitude, end_latitude, EarthRad = EarthRad, Warn = FALSE)
	
	# Get starting new position (that will be overwritten until the collision point is found):
	new_position <- EndPoint(start_longitude, start_latitude, bearing, distance_modifier * distance, EarthRad = EarthRad)[c("long", "lat")]
	
	# Whilst the collision point has not been found:
	while(all.equal(as.vector(GreatCircleDistanceFromLongLat(continent_longitude, continent_latitude, new_position$long, new_position$lat, EarthRad = EarthRad, Warn = FALSE)), continent_radius) != TRUE) {
		
		# Record current distance of new position from centre of continent:
		current_distance_from_centre <- as.vector(GreatCircleDistanceFromLongLat(continent_longitude, continent_latitude, new_position$long, new_position$lat, EarthRad = EarthRad, Warn = FALSE))
		
		# Do we need to increase the distance modifier value?:
		if(current_distance_from_centre < continent_radius) {
			
			# Create potential better value by adding step size to modifier:
			limit <- distance_modifier + stepsize
			
		# Or do we need to decrease the degree modifier value?:
		} else {
			
			# Create potential better value by subtracting step size from modifier:
			limit <- distance_modifier - stepsize
			
		}
		
		# Store the new new position (so we can ask if this is better):
		new_new_position <- EndPoint(start_longitude, start_latitude, bearing, limit * distance, EarthRad = EarthRad)[c("long", "lat")]
		
		# Get new distance from continent centre based on potential better modifer:
		new_distance_from_centre <- as.vector(GreatCircleDistanceFromLongLat(continent_longitude, continent_latitude, new_new_position$long, new_new_position$lat, EarthRad = EarthRad, Warn = FALSE))
		
		# If new distance is closer to continent radius:
		if(abs(current_distance_from_centre - continent_radius) > abs(new_distance_from_centre - continent_radius)) {
			
			# Update new position:
			new_position <- new_new_position
			
			# Update distance modifier itself:
			distance_modifier <- limit
			
		# If current distance is stil our best estimate:
		} else {
			
			# Shrink the step size so we can hone in closer:
			stepsize <- stepsize * 0.1
			
		}
		
	}
	
	# Isolate collision point:
	new_position <- as.vector(unlist(new_position))
	
	# Find remaining (unspent) distance of dispersal step:
	remaining_distance <- as.vector(distance - GreatCircleDistanceFromLongLat(start_longitude, start_latitude, new_position[1], new_position[2], EarthRad = EarthRad, Warn = FALSE))
	
	# Compile output as list:
	output <- list(new_position[1], new_position[2], remaining_distance)
	
	# Add names to output:
	names(output) <- c("collision_longitude", "collision_latitude", "unspent_distance")
	
	# Return output:
	return(output)
	
}
