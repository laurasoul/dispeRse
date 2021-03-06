#' Finds the spherical angle for an equilateral triangle
#'
#' Returns the spherical angle in degrees for an equilateral triangle of known side length
#'
#' @param side_length The great circle distance of a side of the triangle in kilometres.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return Spherical angle in degrees.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' SphericalAngleForEquilateralTriangleFromGreatCircleSideLength(1000)

# Get the spherical angle for an equilateral triangle where its sides (as great circle distances) are known:
SphericalAngleForEquilateralTriangleFromGreatCircleSideLength <- function(side_length, EarthRad = 6367.4447) {
	
	# Error check for large triangles:
	if(side_length > (0.5 * pi * EarthRad)) stop("ERROR: This function will not work for spherical equilateral triangles that are greater in area than 1/8th the sphere's surface.")
	
	# Get y (which will be double the opposite for the triangle we are trying to solve):
	y <- ChordLengthFromTheta(ThetaFromGreatCircleDistance(side_length))
	
	# Get z which will be double the hypotenuse for the triangle we are trying to solve):
	z <- ChordLengthFromTheta(ThetaFromGreatCircleDistance(side_length * 2))
	
	# Solve the traingle for theta (the spherical angle) and convert to degrees:
	spherical_angle <- asin((y / 2) / (z / 2)) / (pi / 180) * 2
	
	# Return the spherical angle:
	return(spherical_angle)
	
}
