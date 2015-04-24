#' Split apart joined continents
#'
#' Return two separate clumps of continents after a separation has occurred.
#'
#' @param min_separation The minimum separation in kilometres between continents.
#' @param longitudes Decimalised longitudes of the continents forming the clump in \code{continent_numbers} order.
#' @param latitudes Decimalised latitudes of the continents forming the clump in \code{continent_numbers} order.
#' @param continent_numbers A character vector of continent numbers.
#' @param protected_links A two-column matrix of protected links (those that cannot be severed).
#' @param EarthRad Earth radius in kilometres.
#' @return Vector of two continental clumps left after separation.
#' @details Nothing yet.
#' @examples
#' # Nothing yet

ContinentSplitter <- function(min_separation, longitudes, latitudes, continent_numbers, protected_links = matrix(nrow=0, ncol=2), EarthRad = 6367.4447) {

	# Check min separation is not greater than 1/2 circumference of planet:
	if(min_separation > (pi * EarthRad)) stop("ERROR: Cannot have minimum separation larger than half the circumference of the planet!")
	
	# Make sure continent numbers are stored as character vector for matching:
	continent_numbers <- as.character(continent_numbers)
	
	# Get all intercontinental links in cluster:
	unprotected_links <- intercontinent_links <- IntercontinentalLinks(min_separation, longitudes, latitudes, EarthRad = EarthRad)
	
	# Add row and column names:
	colnames(unprotected_links) <- colnames(intercontinent_links) <- rownames(unprotected_links) <- rownames(intercontinent_links) <- continent_numbers
	
	# If there are protected links:
	if(nrow(protected_links) > 0) {
		
		# For each protected link:
		for(i in 1:nrow(protected_links)) {
		
			# Update unprotected links by setting protected link state to zero:
			unprotected_links[as.character(protected_links[i, 1]), as.character(protected_links[i, 2])] <- unprotected_links[as.character(protected_links[i, 2]), as.character(protected_links[i, 1])] <- 0
			
		}
		
	}
	
	# If there are unprotected links (and can continue with separation):
	if(sum(unprotected_links[lower.tri(unprotected_links)]) > 0) {
	
		# Pick a random link from the matrix:
		link_position <- sample(which(intercontinent_links == 1))[1]
		
		# Set starting link (from continent to continent):
		start_link <- sort(c(link_position %% nrow(intercontinent_links), ceiling(link_position / nrow(intercontinent_links))))
		
		# Whilst randomly chosen link is a protected link:
		while(unprotected_links[start_link[1], start_link[2]] == 0) {
			
			# Pick a new random link from the matrix:
			link_position <- sample(which(intercontinent_links == 1))[1]
			
			# Update start link:
			start_link <- sort(c(link_position %% nrow(intercontinent_links), ceiling(link_position / nrow(intercontinent_links))))
			
		}
		
		# Define start links by continent numbers:
		start_link <- colnames(intercontinent_links)[start_link]
		
		# Get starting longitude for randomly selected link:
		random_link_start_longitude <- longitudes[match(start_link[1], colnames(intercontinent_links))]
		
		# Get starting latitude for randomly selected link:
		random_link_start_latitude <- latitudes[match(start_link[1], colnames(intercontinent_links))]
		
		# Get ending longitude for randomly selected link:
		random_link_end_longitude <- longitudes[match(start_link[2], colnames(intercontinent_links))]
		
		# Get ending latitude for randomly selected link:
		random_link_end_latitude <- latitudes[match(start_link[2], colnames(intercontinent_links))]
		
		# Get bearing from start of random link to end of random link:
		random_link_start_bearing <- BearingBetweenTwoLongLatPoints(random_link_start_longitude, random_link_start_latitude, random_link_end_longitude, random_link_end_latitude)
		
		# Pick a random distance along the link:
		random_distance <- runif(1, 0, min_separation)
		
		# Define point at which cut will begin:
		cut_point_1 <- EndPoint(random_link_start_longitude, random_link_start_latitude, random_link_start_bearing, random_distance)[c("long", "lat")]
		
		# Get cut bearing (will be used to describe Great Circle made by cut):
		cut_bearing <- runif(1, 0, 360)
		
		# Pick second point on Great Circle to use to describe it:
		cut_point_2 <- EndPoint(cut_point_1$long, cut_point_1$lat, cut_bearing, min_separation)[c("long", "lat")]
		
		# If there are protected links (need to check if cut line intersects them and if so redraw):
		if(nrow(protected_links) > 0) {
			
			# Variable to switch to true if intersection is found:
			intersection_occurs <- FALSE
			
			# For each protected link:
			for(i in 1:nrow(protected_links)) {
				
				# Find if there are any intersections:
				intersections <- ArcIntersection(longitudes[match(as.character(protected_links[i, ])[1], colnames(intercontinent_links))], latitudes[match(as.character(protected_links[i, ])[1], colnames(intercontinent_links))], longitudes[match(as.character(protected_links[i, ])[2], colnames(intercontinent_links))], latitudes[match(as.character(protected_links[i, ])[2], colnames(intercontinent_links))], cut_point_1$long, cut_point_1$lat, cut_point_2$long, cut_point_2$lat, type = c("arc", "GC"), EarthRad = EarthRad)
			
				# If there are intersectiosn update intersection_occurs:
				if(nrow(intersections) > 0) intersection_occurs <- TRUE
				
			}
			
			# Set up counter (to be used rto warn user if loop gets stuck):
			counter <- 1
			
			# Whilst there is an intersection between the cut line and a protected link:
			while(intersection_occurs) {
			
				# Get new cut bearing (will be used to describe Great Circle made by cut):
				cut_bearing <- runif(1, 0, 360)
				
				# Pick new second point on Great Circle to use to describe it:
				cut_point_2 <- EndPoint(cut_point_1$long, cut_point_1$lat, cut_bearing, min_separation)[c("long", "lat")]
				
				# Overwrite intersection occurs:
				intersection_occurs <- FALSE
				
				# For each protected link:
				for(i in 1:nrow(protected_links)) {
					
					# Find if there are any intersections:
					intersections <- ArcIntersection(longitudes[match(as.character(protected_links[i, ])[1], colnames(intercontinent_links))], latitudes[match(as.character(protected_links[i, ])[1], colnames(intercontinent_links))], longitudes[match(as.character(protected_links[i, ])[2], colnames(intercontinent_links))], latitudes[match(as.character(protected_links[i, ])[2], colnames(intercontinent_links))], cut_point_1$long, cut_point_1$lat, cut_point_2$long, cut_point_2$lat, type = c("arc", "GC"), EarthRad = EarthRad)
					
					# If there are intersectiosn update intersection_occurs:
					if(nrow(intersections) > 0) intersection_occurs <- TRUE
					
				}
				
				# Update counter:
				counter <- counter + 1
				
				# Add stop in case this loop never closes:
				if(counter == 10000) stop("ERROR: Model failed as cannot find clean separation cut line.")
			
			}
				
		}
		
		# Get first pole to the cut line equator:
		cut_line_pole_1 <- EndPoint(cut_point_1$long, cut_point_1$lat, (cut_bearing + 90) %% 360, 0.5 * pi * EarthRad)[c("long", "lat")]
		
		# Get second pole to the cut line equator:
		cut_line_pole_2 <- EndPoint(cut_point_1$long, cut_point_1$lat, (cut_bearing - 90) %% 360, 0.5 * pi * EarthRad)[c("long", "lat")]
		
		# Get Great Circle distances from first pole:
		GC_distances_to_pole_1 <- One2ManyGreatCircleDistance(cut_line_pole_1$long, cut_line_pole_1$lat, longitudes, latitudes, EarthRad = EarthRad)

		# Get Great Circle distances from second pole:
		GC_distances_to_pole_2 <- One2ManyGreatCircleDistance(cut_line_pole_2$long, cut_line_pole_2$lat, longitudes, latitudes, EarthRad = EarthRad)
		
		# Make first clump (those in first hemisphere):
		first_clump <- paste(continent_numbers[which(GC_distances_to_pole_1 < (0.5 * pi * EarthRad))], collapse="&")
		
		# Make second clump (those in second hemisphere):
		second_clump <- paste(continent_numbers[which(GC_distances_to_pole_2 < (0.5 * pi * EarthRad))], collapse="&")
		
		# Make vector of clumps:
		clumps <- c(first_clump, second_clump)
		
	# Case if no unprotected links:
	} else {
		
		# Make single clump from all continents:
		clumps <- paste(continent_numbers, collapse="&")
		
		# Warn user that no separation can be made:
		print("WARNING: Cannot separate continent as all links are protected.")
		
	}
	
	# Output clumps:
	return(clumps)
	
}
