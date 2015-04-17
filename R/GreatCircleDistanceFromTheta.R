# Get great circle distacne from theta in radians and radius in km:
GreatCircleDistanceFromTheta <- function(theta, EarthRad = 6367.4447) {
	
	# Calculate great circle distance from theta in radians:
	theta * EarthRad
	
}
