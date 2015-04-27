#' Point(s) at which two arcs intersect on a sphere
#'
#' Given two arcs returns point(s) at which they intersect on a sphere.
#'
#' @param longitude_1 Decimalised longitude of first point on first arc.
#' @param latitude_1 Decimalised latitude of first point on first arc.
#' @param longitude_2 Decimalised longitude of second point on first arc.
#' @param latitude_2 Decimalised latitude of second point on first arc.
#' @param longitude_3 Decimalised longitude of first point on second arc.
#' @param latitude_3 Decimalised latitude of first point on second arc.
#' @param longitude_4 Decimalised longitude of second point on second arc.
#' @param latitude_4 Decimalised latitude of second point on second arc.
#' @param type Vector to indicate if elements are arcs or Great Circles.
#' @param EarthRad Earth radius in kilometres.
#' @return Matrix of longitude-latitude points at which intersection(s) occur.
#' @details Assumes shortest distance between points describes arc. Can define Great Circle-arc intersection with \code{type} option. 
#'
#' @examples
#' longitude_1 <- runif(1, -180, 180)
#' longitude_2 <- runif(1, -180, 180)
#' longitude_3 <- runif(1, -180, 180)
#' longitude_4 <- runif(1, -180, 180)
#' latitude_1 <- runif(1, -90, 90)
#' latitude_2 <- runif(1, -90, 90)
#' latitude_3 <- runif(1, -90, 90)
#' latitude_4 <- runif(1, -90, 90)
#' ArcIntersection(longitude_1, latitude_1, longitude_2, latitude_2,
#'   longitude_3, latitude_3, longitude_4, latitude_4, EarthRad = 6367.4447)

ArcIntersection <- function(longitude_1, latitude_1, longitude_2, latitude_2, longitude_3, latitude_3, longitude_4, latitude_4, type = c("arc", "arc"), EarthRad = 6367.4447) {

	# Check type is of correct length:
	if(length(type) != 2) stop("ERROR: Type must be of length two.")
	
	# Check elements of type are of the correct type:
	if(sum((type == "GC") + (type == "arc")) != 2) stop("ERROR: Elements of type must be either \"arc\" or \"GC\" only.")
	
	# Get intersection points of the two Great Circles described by the arcs:
	res <- GreatCircleIntersection(longitude_1, latitude_1, longitude_2, latitude_2, longitude_3, latitude_3, longitude_4, latitude_4, EarthRad = EarthRad)

	# Only proceed if at least one element is an arc (otherwise res already shows GC-GC intersection pointsand no changes needed):
	if(length(sort(match("arc", type)))) {
	
		# Calculate distance from first point on first arc to first intersection point:
		distance_to_first_intersection_1 <- GreatCircleDistanceFromLongLat(longitude_1, latitude_1, res[, "lon1"], res[, "lat1"], EarthRad = EarthRad)
		
		# Calculate distance from first point on first arc to second intersection point:
		distance_to_second_intersection_1 <- GreatCircleDistanceFromLongLat(longitude_1, latitude_1, res[, "lon2"], res[, "lat2"], EarthRad = EarthRad)
		
		# Calculate distance from first point on second arc to first intersection point:
		distance_to_first_intersection_2 <- GreatCircleDistanceFromLongLat(longitude_3, latitude_3, res[, "lon1"], res[, "lat1"], EarthRad = EarthRad)
		
		# Calculate distance from first point on second arc to second intersection point:
		distance_to_second_intersection_2 <- GreatCircleDistanceFromLongLat(longitude_3, latitude_3, res[, "lon2"], res[, "lat2"], EarthRad = EarthRad)

		# If first element is an arc:
		if(type[1] == "arc") {
			
			# Set distance as smallest Great Circle distance between the two points:
			distance_to_other_point_1 <- GreatCircleDistanceFromLongLat(longitude_1, latitude_1, longitude_2, latitude_2, EarthRad = EarthRad)
		
		# If first element is a Great Circle:
		} else {
		
			# Set distance as full circumference of planet:
			distance_to_other_point_1 <- 2 * pi * EarthRad
			
		}

		# If second element is an arc:
		if(type[2] == "arc") {
			
			# Set distance as smallest Great Circle distance between the two points:
			distance_to_other_point_2 <- GreatCircleDistanceFromLongLat(longitude_3, latitude_3, longitude_4, latitude_4, EarthRad = EarthRad)
		
		# If second element is a Great Circle:
		} else {
	
			# Set distance as full circumference of planet:
			distance_to_other_point_2 <- 2 * pi * EarthRad
			
		}
			
		# Create empty delete rows vector (will store intersection points outside of arcs):
		delete.rows <- vector(mode="numeric")
		
		# Make result into matrix:
		res <- matrix(res, ncol=2, byrow=TRUE, dimnames=list(c(), c("Longitude", "Latitude")))
		
		# Case if first point is not an actual intersection:
		if(!(distance_to_other_point_1 >= distance_to_first_intersection_1 && distance_to_other_point_2 >= distance_to_first_intersection_2)) delete.rows <- c(delete.rows, 1)
		
		# Case if second point is not an actual intersection:
		if(!(distance_to_other_point_1 >= distance_to_second_intersection_1 && distance_to_other_point_2 >= distance_to_second_intersection_2)) delete.rows <- c(delete.rows, 2)
		
		# Remove intersections outside of arcs:
		res <- res[-delete.rows, ]
		
		# Ensure res is still a matrix:
		if(!is.matrix(res)) res <- matrix(res, ncol=2, byrow=TRUE, dimnames=list(c(), c("Longitude", "Latitude")))
	
	# Case if a Great Circle-Great Circle intersection:
	} else {
		
		# Make result into matrix:
		res <- matrix(res, ncol=2, byrow=TRUE, dimnames=list(c(), c("Longitude", "Latitude")))
		
	}
		
	# Return matrix of zero, one, or two rows (depending on number of intersection points):
	return(res)
	
}
