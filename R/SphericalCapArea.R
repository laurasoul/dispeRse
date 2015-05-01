#' Spherical Cap Area
#'
#' Returns the surface area of a cap (circle) on a sphere
#'
#' @param radius Radius of the cap in kilometres.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return The surface area of the cap in square-kilometres.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' SphericalCapArea(2000)

SphericalCapArea <- function(radius, EarthRad = 6367.4447) {
	
	# Get circumference of the sphere:
	circumference <-  2 * pi * EarthRad
	
	# Get the degrees of latitude to which this corresponds:
	degrees_of_latitude <- radius / circumference * 360
	
	# Use Pythagoras to get the radius of the sphere minus the height of the cap defined by the circle:
	radius_minus_height <- cos(degrees_of_latitude * (pi / 180)) * EarthRad
	
	# Now get the height of the cap:
	height <- EarthRad - radius_minus_height
	
	# Get the cap area:
	cap_area <- 2 * pi * EarthRad * height
	
	# Return the cap area:
	return(cap_area)
	
}
