# Get chord length from theta in radians and radius in km:
ChordLengthFromTheta <- function(theta, EarthRad = 6367.4447) {
	
	# Get chord length from theta (in radians):
	sin(theta / 2) * EarthRad * 2
	
}
