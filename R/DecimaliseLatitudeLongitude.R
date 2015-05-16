#' Decimalises a longitude or latitude
#'
#' Decimalises a longitude or latitude specified in degrees, minutes, seconds, and decimals
#'
#' @param degrees An integer specifying the longitude or latitiude in degrees.
#' @param minutes An integer specifying the longitude or latitiude minutes.
#' @param seconds An integer specifying the longitude or latitiude seconds.
#' @param decimals An integer specifying the longitude or latitiude decimals.
#' @param direction The direction (N = North, E = East, S = South, W = West).
#' @return A decimalised latitude or longitude.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' DecimaliseLatitudeLongitude(degrees = 179, minutes = 12,
#'   seconds = 13, decimals = 67, direction = "W")

DecimaliseLatitudeLongitude <- function(degrees, minutes = 0, seconds = 0, decimals = 0, direction) {

	# Check direction input is correct:
	if(direction != "N" && direction != "E" && direction != "S" && direction != "W") stop("ERROR: Direction must be one of \"N\", \"S\", \"E\", or \"W\".")
	
	# Case if input is a latitude:
	if(direction == "N" || direction == "S") {
		
		# Check input is between 0 and 90:
		if(degrees > 90 || degrees < 0) stop("ERROR: Latitude must be between 0 and 90 degrees.")
		
		# Check input does not exceed 90:
		if(degrees == 90 && (minutes + seconds + decimals) > 0) stop("ERROR: Latitude must be between 0 and 90 degrees.")
		
	# Case if input is a longitude:
	} else {
		
		# Check input is between 0 and 180:
		if(degrees > 180 || degrees < 0) stop("ERROR: Longitude must be between 0 and 180 degrees.")

		# Check input does not exceed 180:
		if(degrees == 180 && (minutes + seconds + decimals) > 0) stop("ERROR: Longitude must be between 0 and 180 degrees.")
		
	}
	
	# Check minutes are in correct format:
	if(minutes > 60 || minutes < 0) stop("ERROR: Minutes must be between 0 and 60.")

	# Check seconds are in correct format:
	if(seconds > 60 || seconds < 0) stop("ERROR: Seconds must be between 0 and 60.")
	
	# Check decimals are in correct format:
	if(decimals > 100 || decimals < 0) stop("ERROR: Decimals must be between 0 and 100.")
	
	# Decimalise minutes:
	minutes <- minutes / 60
	
	# Decimalise seconds:
	seconds <- seconds / 6000
	
	# Decimalise decimals:
	decimals <- decimals / 1000000

	# Create latitude output:
	if(direction == "E" || direction == "N") out <- degrees + minutes + seconds + decimals
	
	# Create longitude output:
	if(direction == "W" || direction == "S") out <- (degrees + minutes + seconds + decimals) * -1
	
	# Return output:
	return(out)

}
