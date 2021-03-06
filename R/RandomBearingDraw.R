#' Make random draws for bearings
#'
#' Make random draws for bearings allowing bias towards North-South or East-West
#' @param n The number of randomly drawn bearings required.
#' @param shape_parameter The shape parameter used as both values in a beta distribution.
#' @return A vector of bearing(s) of length n.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
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
