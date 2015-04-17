# Get great circle distacne from theta in radians and radius in km:
GetGreatCircleDistanceFromTheta <- function(theta, EarthRad = 6367.4447) {
	
	# Calculate great circle distance from theta in radians:
	theta * EarthRad
	
}
