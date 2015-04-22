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
#'   continent_speed_mean = 500, continent_speed_sd = 250, EarthRad = 6367.4447,
#'   polar = FALSE)

# Outputs:
# - N separated continents (distance matrix with values less than 2 radii)
# - Long-lat of each circle centre
# - Bearings after each step change
# - Total land area all circles - minus

EverythingFunction <- function(N_steps = 1000, N_continents = 7, radius = 2000, start_configuration = "supercontinent", squishiness = 0.25, stickiness = 0.95, continent_speed_mean = 500, continent_speed_sd = 200, EarthRad = 6367.4447, polar = FALSE) {

# Need more top-level conditionals, e.g. N steps must be a positive integer, speed mean and sd must also be positive
# Others may be cuaght by subfunctions so no need to repeat
	
	# Start by picking continent start points:
	continent_starting_points <- StartingPoints(N_continents = N_continents, radius = radius, start_configuration = start_configuration, squishiness = squishiness, EarthRad = EarthRad, polar = polar)
	
# Are any continents joined at start?:

	# Get minimum_continental separation:
	min_separation <- (1 - squishiness) * radius * 2

	# Get list of separate continents:
	separate_continents <- HowManySeparateContinents(min_separation, continent_starting_points[, "Longitude"], continent_starting_points[, "Latitude"])
	
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
		euler_GC_distances <- GreatCircleDistanceMatrix(rbind(c(euler_pole_longitudes[i], euler_pole_latitudes[i]), continent_starting_points[, c("Longitude", "Latitude")])[, "Longitude"], rbind(c(euler_pole_longitudes[i], euler_pole_latitudes[i]), continent_starting_points[, c("Longitude", "Latitude")])[, "Latitude"])[2:(N_continents + 1), 1]
		
		# Find GC distance to furthest continent (closest to euler pole "equator") in cluster (as speed will be assigned based on this):
		furthest_continent_GC_distance <- euler_GC_distances[which.min(abs(euler_GC_distances - rep(0.5 * pi * EarthRad, length(euler_GC_distances))))]
		
		# Randomly draw a continent speed:
		continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
		# If a negative or zero speed is picked then redraw:
		while(continent_speed <= 0) continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
		# Set degree change per step (effectively the speed):
		degrees_per_step[i] <- continent_speed / (2 * pi * furthest_continent_GC_distance) * 360
		
	}
	
# Now need to move them!
	
	
	
# When rotating around Euler pole could theoretically pick clockwise or anticlockwise, but as we are allowing poles to be on either side of planet this takes care of that for us!
# Number continents in plots
	
}
