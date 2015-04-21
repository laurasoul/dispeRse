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
#' @return A magic table of awesomeness.
#' @details Nothing yet.
#'
#' @examples
#' EverythingFunction <- function(N_steps = 1000, N_continents = 7, radius = 2000,
#'   start_configuration = "supercontinent", squishiness = 0.25, stickiness = 0.95,
#'   continent_speed_mean = 500, continent_speed_sd = 250, EarthRad = 6367.4447)

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

EverythingFunction <- function(N_steps = 1000, N_continents = 7, radius = 2000, start_configuration = "supercontinent", squishiness = 0.25, stickiness = 0.95, continent_speed_mean = 500, continent_speed_sd = 250, EarthRad = 6367.4447, polar = FALSE) {

	# Start by picking continent start points:
	continent_starting_points <- StartingPoints(N_continents = N_continents, radius = radius, start_configuration = start_configuration, squishiness = squishiness, EarthRad = EarthRad, polar = polar)
	
	
# Are any continents joined at start?:

	# Get minimum_continental separation:
	min_separation <- (1 - squishiness) * radius * 2

	HowManySeparateContinents(min_separation, continent_starting_points[, "Longitude"], continent_starting_points[, "Latitude"])
	
}
