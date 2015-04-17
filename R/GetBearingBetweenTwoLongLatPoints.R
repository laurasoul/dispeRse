# Get bearing (with 0 or 360 being North, 90 being East etc.) from one lat-long coordinate to another:
GetBearingBetweenTwoLongLatPoints <- function(start_long, start_lat, end_long, end_lat) {
	
# Check that lat and long are not both identical (where no bearing is possible)
	
	# Convert starting longitude to radians:
	start_long <- start_long * (pi / 180)
	
	# Convert starting latitude to radians:
	start_lat <- start_lat * (pi / 180)
	
	# Convert ending longitude to radians:
	end_long <- end_long * (pi / 180)
	
	# Convert ending latitude to radians:
	end_lat <- end_lat * (pi / 180)
	
	# Get bearing in radians:
	bearing <- atan2(sin(end_long - start_long) * cos(end_lat), cos(start_lat) * sin(end_lat) - sin(start_lat) * cos(end_lat) * cos(end_long - start_long))
	
	# Convert bearing to degrees:
	bearing <- bearing / (pi / 180)
	
	# Place bearing on 0 to 360 scale, modulo 360:
	bearing <- (bearing + 360) %% 360
	
	# Return bearing:
	return(bearing)
	
}
