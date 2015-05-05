#' Sweepstakes dispersal
#'
#' If sweeptakes dispersal occurs picks a landing spot.
#'
#' @param start_longitude Decimalised longitude of organism's dispersal starting point.
#' @param start_latitude Decimalised latitude of organism's dispersal starting point.
#' @param end_longitude Decimalised longitude of organism's dispersal ending point.
#' @param end_latitude Decimalised latitude of organism's dispersal ending point.
#' @param start_continent Number of starting continent.
#' @param touching_continents Character vector of continental clusters separated by ampersands.
#' @param continent_centres A two-column (longitude, latitude) of the centres of each continent in order (1 at top to N at bottom).
#' @param continent_radius The radius of the continent in kilometres.
#' @param EarthRad Earth radius in kilometres.
#' @return The new position of the organism following a sweepstakes dispersal.
#' @author Laura C. Soul \email{lauracsoul@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' #MagicPortal()

# Maybe limit maximum dispersal distance? Currently this is effectively halfway round the planet.

MagicPortal <- function(start_longitude, start_latitude, end_longitude, end_latitude, start_continent, touching_continents, continent_centres, continent_radius, EarthRad = 6367.4447) {
	
	# Create empty variables to store landing spot and landing continent:
	landing_spot <- landing_continent <- NA
	
	# Case if there are no other continents (i.e., there is a single supercontinent):
	if(length(touching_continents) == 1) {
		
		# Cannot disperse to new continent, but can disperse along coast of current continent:
		next
		
	# Case if there is at least one other continental cluster:
	} else {
		
		# Get vector of other continental clusters:
		other_continents <- touching_continents[-WhichSupercontinent(start_continent, touching_continents)]
		
		# Find collision point:
		collision_point <- DivingBoard(start_longitude, start_latitude, end_longitude, end_latitude, continent_centres[as.numeric(start_continent), 1], continent_centres[as.numeric(start_continent), 2], continent_radius, EarthRad = EarthRad)
		
		# Find bearing from continent centre to collision point:
		bearing_to_collision <- BearingBetweenTwoLongLatPoints(continent_centres[as.numeric(start_continent), 1], continent_centres[as.numeric(start_continent), 2], collision_point$collision_longitude, collision_point$collision_latitude)
		
		# Find bearing from collision point out to sea:
		bearing_from_collision <- EndPoint(continent_centres[as.numeric(start_continent), 1], continent_centres[as.numeric(start_continent), 2], bearing_to_collision, continent_radius)$bearing
		
		# Create empty vectors to store bearings and distances to other each continent not in the current cluster:
		other_continent_bearings <- other_continent_distances <- vector(mode="numeric")
		
		# For each such continent:
		for(i in as.numeric(unlist(strsplit(other_continents, "&")))) {
			
			# Find and store the shortest Great Circle distance:
			other_continent_distances <- c(other_continent_distances, GreatCircleDistanceFromLongLat(collision_point$collision_longitude, collision_point$collision_latitude, continent_centres[i, 1], continent_centres[i, 2]))
			
			# Find and store the bearing along the shortest Great Circle distance:
			other_continent_bearings <- c(other_continent_bearings, BearingBetweenTwoLongLatPoints(collision_point$collision_longitude, collision_point$collision_latitude, continent_centres[i, 1], continent_centres[i, 2]))
			
			# Add number of continent for later retrieval:
			names(other_continent_bearings)[length(other_continent_bearings)] <- names(other_continent_distances)[length(other_continent_distances)] <- i

		}

		# Rotate bearings so that the bearing from collision is at 180:
		rotated_bearings <- (other_continent_bearings + (180 - bearing_from_collision)) %% 360
		
		# Get all continents found within plus/minus 90-degrees of collision point:
		potential_landing_spots <- names(other_continent_bearings)[which((as.numeric(rotated_bearings <= 270) + as.numeric(rotated_bearings >= 90)) == 2)]
		
		# Case if there is a potential landing spot:
		if(length(potential_landing_spots) > 0) {
			
			landing_continent <- names(which.min(other_continent_distances[potential_landing_spots]))
			
			landing_spot <- EndPoint(continent_centres[as.numeric(landing_continent), 1], continent_centres[as.numeric(landing_continent), 2], BearingBetweenTwoLongLatPoints(continent_centres[as.numeric(landing_continent), 1], continent_centres[as.numeric(landing_continent), 2], collision_point$collision_longitude, collision_point$collision_latitude), continent_radius)[c("long", "lat")]
			
		# Case if there are no potential landing spots on other continent clusters:
		} else {
			
			# Cannot disperse to new continent, but can disperse along coast of current continent:
			next
			
		}
		
	}
	
	# Case if dispersing onto current continental cluster (no viable landing spots elsewhere):
	if(is.na(landing_spot)) {
		
		# Find current cluster:
		current_cluster <- touching_continents[WhichSupercontinent(start_continent, touching_continents)]
		
		# Isolate component continents:
		current_cluster <- as.numeric(unlist(strsplit(current_cluster, "&")))
		
		# Randomly draw a landing continent:
		landing_continent <- sample(x = rep(current_cluster, 2), size = 1)
		
		# Randomly draw a landing spot (on the coast):
		landing_spot <- EndPoint(continent_centres[landing_continent, 1], continent_centres[landing_continent, 2], runif(1, min = 0, max = 360), continent_radius)
		
		# Calculate distances from landing spot to continent centre(s) (to be used to establish if spot is on the coast or not):
		distances <- One2ManyGreatCircleDistance(landing_spot$long, landing_spot$lat, continent_centres[current_cluster, 1], continent_centres[current_cluster, 2], EarthRad = EarthRad)
		
		# If current landing spot is not on the coast:
		while((as.numeric(apply(matrix(distances, nrow=1), 2, all.equal, current = continent_radius) != TRUE) + as.numeric(distances < continent_radius)) == 2) {
			
			# Redraw landing continent:
			landing_continent <- sample(x = rep(current_cluster, 2), size = 1)
			
			# Redraw landing spot:
			landing_spot <- EndPoint(continent_centres[landing_continent, 1], continent_centres[landing_continent, 2], runif(1, min = 0, max = 360), continent_radius)
			
			# Recalculate distances to continent centre(s):
			distances <- One2ManyGreatCircleDistance(landing_spot$long, landing_spot$lat, continent_centres[current_cluster, 1], continent_centres[current_cluster, 2], EarthRad = EarthRad)

		}
		
	}
	
	# Compile output:
	output <- list(as.numeric(landing_continent), as.vector(unlist(landing_spot)[1]), as.vector(unlist(landing_spot)[2]))
	
	# Add names to output:
	names(output) <- c("continent", "longitude", "latitude")
	
	# Return output:
	return(output)
	
}
