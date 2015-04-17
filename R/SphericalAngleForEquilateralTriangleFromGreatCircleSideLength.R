#' Finds the spherical angle for an equilateral triangle
#'
#' Returns the spherical angle in degrees for an equilateral triangle of known side length
#'
#' @param side_length The great circle distance of a side of the triangle in kilometres.
#' @param EarthRad Radius of the Earth in kilometres.
#' @return Spherical angle in degrees.
#' @details Nothing yet.
#' @examples
#' SphericalAngleForEquilateralTriangleFromGreatCircleSideLength(1000)

# Get the spherical angle for an equilateral triangle where its sides (as great circle distances) are known:
SphericalAngleForEquilateralTriangleFromGreatCircleSideLength <- function(side_length, EarthRad = 6367.4447) {
	
# Check this works for triangles larger than equivalent of 90-degrees of latitude on each side!
	
	# Get y (which will be double the opposite for the triangle we are trying to solve):
	y <- GetChordLengthFromTheta(GetThetaFromGreatCircleDistance(side_length))
	
	# Get z which will be double the hypotenuse for the triangle we are trying to solve):
	z <- GetChordLengthFromTheta(GetThetaFromGreatCircleDistance(side_length * 2))
	
	# Solve the traingle for theta (the spherical angle) and convert to degrees:
	spherical_angle <- asin((y / 2) / (z / 2)) / (pi / 180) * 2
	
	# Return the spherical angle:
	return(spherical_angle)
	
}
