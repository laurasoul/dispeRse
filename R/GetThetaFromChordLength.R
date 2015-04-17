# Get theta in radians from chord length and radius in km:
GetThetaFromChordLength <- function(chord_length, EarthRad = 6367.4447) {
	
	# Get theta in radians from chord length and radius in km:
	asin((chord_length * 0.5) / EarthRad) * 2
	
}
