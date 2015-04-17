# Get theta in radians from great circle distance:
GetThetaFromGreatCircleDistance <- function(great_circle_distance, EarthRad = 6367.4447) {
	
	# Get theta in radians from great circle distance:
	great_circle_distance / EarthRad
	
}
