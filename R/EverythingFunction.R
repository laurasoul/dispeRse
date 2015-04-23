#' Global function to simulate evolution on a sphere
#'
#' Global function to generate simulated continents and clades on a sphere
#'
#' @param N_steps The number of time steps in the simulation.
#' @param N_continents The (maximum) number of individual continents.
#' @param radius The radius of each circular continent.
#' @param start_configuration One of "random separate", "random overlap", "supercontinent", or "max separate".
#' @param squishiness A value from 0 (continents can never overlap) to 1 (continents can overlap completely).
#' @param stickiness Probability of two conjoined continents remaining conjoined in the next time step.
#' @param continent_speed_mean Mean continent speed (kilometers per time step).
#' @param continent_speed_sd Standard deviation of continent speed (kilometers per time step).
#' @param EarthRad Radius of the Earth in kilometres.
#' @param polar Whether continents are allowed to exist exactly at the poles.
#' @return A magic table of awesomeness.
#' @details Nothing yet.
#'
#' @examples
#' EverythingFunction <- function(N_steps = 1000, N_continents = 7, radius = 2000,
#'   start_configuration = "supercontinent", squishiness = 0.25, stickiness = 0.95,
#'   continent_speed_mean = 5, continent_speed_sd = 2, EarthRad = 6367.4447,
#'   polar = FALSE)

# Outputs:
# - N separated continents (distance matrix with values less than 2 radii)
# - Long-lat of each circle centre
# - Bearings after each step change
# - Total land area all circles - minus

EverythingFunction <- function(N_steps = 1000, N_continents = 7, radius = 2000, start_configuration = "supercontinent", squishiness = 0.25, stickiness = 0.95, continent_speed_mean = 5, continent_speed_sd = 2, EarthRad = 6367.4447, polar = FALSE) {

	#Subfunction to find out which supercontinent a circle belongs to
	which_sprcont <- function(cont, sprconts){
		cont <- as.character(cont)
		result <- which(unlist(lapply(lapply(lapply(strsplit(sprconts, "&"), match, cont), sort), length)) == 1)
		return(result)
	}
# Need more top-level conditionals, e.g. N steps must be a positive integer, speed mean and sd must also be positive
# Others may be cuaght by subfunctions so no need to repeat
	
	# Start by picking continent start points:
	continent_starting_points <- StartingPoints(N_continents = N_continents, radius = radius, start_configuration = start_configuration, squishiness = squishiness, EarthRad = EarthRad, polar = polar)
	
# Are any continents joined at start?:

	# Get minimum_continental separation:
	min_separation <- (1 - squishiness) * radius * 2

	# Get list of separate continents:
	separate_continents <- HowManySeparateContinents(min_separation, continent_starting_points[, "Longitude"], continent_starting_points[, "Latitude"])

	# Get list of touching continents (to be used later for whether dispersal is allowable or not):
	touching_continents <- HowManySeparateContinents((radius * 2), continent_starting_points[, "Longitude"], continent_starting_points[, "Latitude"])
	
# Assign Euler poles to each separate continent:
	
	# Randomly draw longitudes for each separated continent:
	euler_pole_longitudes <- runif(length(separate_continents), -180, 180)
	
	# Randomly draw latitudes for each separated continent:
	euler_pole_latitudes <- runif(length(separate_continents), -90, 90)
	
	# Ensure no latitude is directly at a pole (North or South) by redrawing if any are found:
	if((sum(euler_pole_latitudes == 90) + sum(euler_pole_latitudes == -90)) > 0) euler_pole_latitudes <- runif(length(separate_continents), -90, 90)

# Assign speeds to each separate continent:
	
	# Create empty vector to store degrees per step (effectively the speed of movement) for each continent:
	degrees_per_step <- vector(mode="numeric")
	
	# For each separate continent:
	for(i in 1:length(separate_continents)) {
		
		# Get Greate Circle distances from Euler pole to each continent centre:
		euler_GC_distances <- One2ManyGreatCircleDistance(euler_pole_longitudes[i], euler_pole_latitudes[i], continent_starting_points[as.numeric(unlist(strsplit(separate_continents[i], "&"))), "Longitude"], continent_starting_points[as.numeric(unlist(strsplit(separate_continents[i], "&"))), "Latitude"])
		
		# Find GC distance to furthest continent (closest to euler pole "equator") in cluster (as speed will be assigned based on this):
		furthest_continent_GC_distance <- max(abs(euler_GC_distances - rep(0.5 * pi * EarthRad)))
		
		# Randomly draw a continent speed:
		continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
		# If a negative or zero speed is picked then redraw:
		while(continent_speed <= 0) continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
		# Set degree change per step (effectively the speed):
		degrees_per_step[i] <- continent_speed / (2 * pi * furthest_continent_GC_distance) * 360
		
	}
	
# Now need to move them!
	#List to store which circles are in each supercontinent (new element added only when it changes)
	linked <- list()
	linked[[1]] <- separate_continents

	#List to store which circles are in each supercontinent (new element added only when it changes)
	touching <- list()
	touching[[1]] <- touching_continents

	#Array to store the positions of every continent at each time step
	position <- array(NA, c(N_continents, N_steps + 1, 2), c("continent", "timestep", "coordinate"))
	position[,1,1] <- continent_starting_points[,"Longitude"]
	position[,1,2] <- continent_starting_points[,"Latitude"]

	#for loops to move everything
	for (t in 2:(N_steps + 1)) {
		for (k in 1:N_continents) {
			#find current longlat of the circle k
			start_long <- position[k, t-1, 1]
			start_lat <- position[k, t-1, 2]

			#Identify which supercontinent, and therefore which element of the euler pole and speed vectors, the circle k belongs to
			where <- which_sprcont(k, tail(linked, n=1)[[1]])

			#Find distance of circle from pole
			distance <- GreatCircleDistanceFromLongLat(long1=start_long,lat1=start_lat, long2=euler_pole_longitudes[where], lat2=euler_pole_latitudes[where])
			
			#Find bearing of circle from pole
			init_bearing <- BearingBetweenTwoLongLatPoints(euler_pole_longitudes[where], euler_pole_latitudes[where], start_long, start_lat)

			#Find the new bearing of circle from pole according to the speed specified
			new_bearing <- (init_bearing + degrees_per_step[where]) %% 360

			#Find the new location of the circle
			new_loc <- EndPoint(euler_pole_longitudes[where], euler_pole_latitudes[where], new_bearing, distance)

			#Add the new loction to the position matrix
			position[k,t,1] <- new_loc$long
			position[k,t,2] <- new_loc$lat
		}
	}
# When rotating around Euler pole could theoretically pick clockwise or anticlockwise, but as we are allowing poles to be on either side of planet this takes care of that for us!
# Number continents in plots
	
}
