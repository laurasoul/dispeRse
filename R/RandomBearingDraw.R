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
#' # Set seed for example:
#' set.seed(21)
#'
#' # Unbiased random draw (N-S and E-W draws equally likely):
#' hist(RandomBearingDraw(n = 10000, shape_parameter = 1), xlim=c(0, 360),
#'   breaks=36, xlab="Bearing (degrees)")
#'
#' # Random draw biased towards E-W draws:
#' hist(RandomBearingDraw(n = 10000, shape_parameter = 10), xlim=c(0, 360),
#'   breaks=36,xlab="Bearing (degrees)")
#'
#' # Random draw biased towards N-S:
#' hist(RandomBearingDraw(n = 10000, shape_parameter = 0.1), xlim=c(0, 360),
#'   breaks=36, xlab="Bearing (degrees)")

RandomBearingDraw <- function(n = 1, shape_parameter = 1) {

	# Make n random draws using shape parameter:
	draws <- (rbeta(n, shape_parameter, shape_parameter) * (sample(c(-180, 180), n, replace=TRUE))) %% 360
	
	# Output n random draw(s):
	return(draws)
	
}