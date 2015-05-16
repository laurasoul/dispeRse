#' Split apart joined continents
#'
#' Return two separate clumps of continents after a separation has occurred.
#' @param min_separation The minimum separation in kilometres between continents.
#' @param longitudes Decimalised longitudes of the continents forming the clump in \code{continent_numbers} order.
#' @param latitudes Decimalised latitudes of the continents forming the clump in \code{continent_numbers} order.
#' @param continent_numbers A character vector of continent numbers.
#' @param protected_links A two-column matrix of protected links (those that cannot be severed).
#' @param EarthRad Earth radius in kilometres.
#' @param Warn Whether or not to print warnings.
#' @return Vector of two continental clumps left after separation.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' # Nothing yet

ContinentSplitter <- function(min_separation, longitudes, latitudes, continent_numbers, protected_links = matrix(nrow = 0, ncol = 2), EarthRad = 6367.4447, Warn = TRUE) {

	# Check min separation is not greater than 1/2 circumference of planet:
	if(min_separation > (pi * EarthRad)) stop("ERROR: Cannot have minimum separation larger than half the circumference of the planet!")
	
	# Make sure continent numbers are stored as character vector for matching:
	continent_numbers <- as.character(continent_numbers)
	
	# Get all intercontinental links in cluster:
	unprotected_links <- intercontinent_links <- IntercontinentalLinks(min_separation, longitudes, latitudes, EarthRad = EarthRad)
	
	# Add row and column names:
	colnames(unprotected_links) <- colnames(intercontinent_links) <- rownames(unprotected_links) <- rownames(intercontinent_links) <- continent_numbers
	
	# Create empty vector to store protected clumps:
	protected_clumps <- vector(mode="character")
	
	# If there are protected links:
	if(nrow(protected_links) > 0) {
		
		# For each protected link:
		for(i in 1:nrow(protected_links)) {
		
			# Update unprotected links by setting protected link state to zero:
			unprotected_links[as.character(protected_links[i, 1]), as.character(protected_links[i, 2])] <- unprotected_links[as.character(protected_links[i, 2]), as.character(protected_links[i, 1])] <- 0
			
		}
		
		# If there are possible additional protected links (by implication as they would lead to impossible positions to find a cut line from):
		if(nrow(protected_links) > 1) {
		
			# For each protected link excep the last:
			for(i in 1:(nrow(protected_links) - 1)) {
				
				# For each next protected link to the last:
				for(j in (i + 1):nrow(protected_links)) {
					
					# Get continents forming the two links:
					continents <- rle(sort(as.vector(protected_links[c(i, j), ])))
					
					# If any continents appear twice (i.e., are links linked?):
					if(any(continents$lengths == 2)) {
						
						# Get ends of linked-links (to check if these form a currently unprotected link that should really be considered protected):
						continents <- as.character(continents$values[continents$lengths == 1])
						
						# If link exists and is unprotected:
						if(unprotected_links[continents[1], continents[2]] == 1) {
							
							# Remove link from unprotected links matrix (i.e., set to zero):
							unprotected_links[continents[1], continents[2]] <- unprotected_links[continents[2], continents[1]] <- 0
							
							# Add new protected links to list:
							protected_links <- rbind(protected_links, sort(as.numeric(c(continents[2], continents[1]))))
							
						}
						
					}
					
				}
				
			}

		}
		
		# For each protected link:
		for(i in 1:nrow(protected_links)) {
			
			# Search for matches of continents involved in link with clumps found so far:
			clump_matches <- unique(c(WhichSupercontinent(protected_links[i, 1], protected_clumps), WhichSupercontinent(protected_links[i, 2], protected_clumps)))
			
			# If no matching clump, then make new clump of continents:
			if(length(clump_matches) == 0) protected_clumps <- c(protected_clumps, paste(sort(protected_links[i, ]), collapse="&"))
			
			# If matching clump found add new continent(s) to it:
			if(length(clump_matches) == 1) protected_clumps[clump_matches] <- paste(sort(unique(c(as.numeric(unlist(strsplit(protected_clumps[clump_matches], "&"))), protected_links[i, ]))), collapse="&")
			
		}
		
	}
	
	# If there are unprotected links (and can continue with separation):
	if(sum(unprotected_links[lower.tri(unprotected_links)]) > 0) {
		
		# If all continents are found in protected clumps (then clumps define the answer and splitting must occur based on these):
		if(length(unlist(strsplit(protected_clumps, "&"))) == length(continent_numbers)) {
			
			# Make protected clumps the output:
			clumps <- protected_clumps
			
			# If only one continent exists in the protected clumps then warn user that no separation has been made:
			if(length(protected_clumps) == 1 && Warn) print("WARNING: Cannot separate continent as protected links prohibit it.")
			
		}
		
		# If only one continent is missing from the protected clump:	
		if((length(unlist(strsplit(protected_clumps, "&"))) + 1) == length(continent_numbers)) {
			
			# Form output clumps from protected clumps and isolated continent:
			clumps <- c(protected_clumps, setdiff(continent_numbers, unlist(strsplit(protected_clumps, "&"))))
			
		}
		
		# Case if two or more continents are missing from protected clumps (actually draw a cut line):
		if((length(unlist(strsplit(protected_clumps, "&"))) + 1) < length(continent_numbers)) {

			# If protected clumps is not empty:
			if(length(protected_clumps) > 0) {

				# Make unprotected links into a matrix of continent-continent links:
				unprotected_links_matrix <- matrix(as.numeric(unlist(strsplit(unique(apply(t(apply(cbind(colnames(unprotected_links)[which(intercontinent_links == 1, arr.ind = TRUE)[, 1]], colnames(unprotected_links)[which(intercontinent_links == 1, arr.ind = TRUE)[, 2]]), 1, sort)), 1, paste, collapse="&")), "&"))), ncol=2, byrow=TRUE)
				
				# Use protected clumps as starting point:
				unprotected_clumps <- protected_clumps
				
				# For each unprotected link:
				for(i in 1:nrow(unprotected_links_matrix)) {
					
					# Search for matches of continents involved in link with clumps found so far:
					clump_matches <- unique(c(WhichSupercontinent(unprotected_links_matrix[i, 1], unprotected_clumps), WhichSupercontinent(unprotected_links_matrix[i, 2], unprotected_clumps)))
					
					# If no matching clump, then make new clump of continents:
					if(length(clump_matches) == 0) unprotected_clumps <- c(unprotected_clumps, paste(sort(unprotected_links_matrix[i, ]), collapse="&"))
					
					# If matching clump found add new continent(s) to it:
					if(length(clump_matches) == 1) unprotected_clumps[clump_matches] <- paste(sort(unique(c(as.numeric(unlist(strsplit(unprotected_clumps[clump_matches], "&"))), unprotected_links_matrix[i, ]))), collapse="&")
					
				}
				
# Remove protected clumps from unprotected clumps
				
				# For each protected clump:
				for(i in 1:length(protected_clumps)) {
					
					# For each continent in a protected clump:
					for(j in unlist(strsplit(protected_clumps[i], "&"))) {
						
						# Which unprotected clump is the continetn found in:
						where <- WhichSupercontinent(j, unprotected_clumps)
						
						# Remove protected clump continent from unprotected clump continent:
						unprotected_clumps[where] <- paste(unlist(strsplit(unprotected_clumps[where], "&"))[-match(j, unlist(strsplit(unprotected_clumps[where], "&")))], collapse="&")
						
					}
					
				}
				
				# Make clumps out of unprotected and protected continental clumps:
				clumps <- c(protected_clumps, unprotected_clumps)[nchar(c(protected_clumps, unprotected_clumps)) > 0]
				
			}
			
			# Case if protected clumps is still empty:
			if(length(protected_clumps) == 0) {
				
				# Pick random starting link (from continent to continent):
				start_link <- which(unprotected_links == 1, arr.ind = TRUE)[sample(1:sum(intercontinent_links == 1))[1],]
				
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
				cut_point_1 <- EndPoint(random_link_start_longitude, random_link_start_latitude, random_link_start_bearing, random_distance, EarthRad = EarthRad)[c("longitude", "latitude")]
				
				# Get cut bearing (will be used to describe Great Circle made by cut):
				cut_bearing <- runif(1, 0, 360)
				
				# Pick second point on Great Circle to use to describe it:
				cut_point_2 <- EndPoint(cut_point_1$longitude, cut_point_1$latitude, cut_bearing, min_separation, EarthRad = EarthRad)[c("longitude", "latitude")]
				
				# If there are protected links (need to check if cut line intersects them and if so redraw):
				if(nrow(protected_links) > 0) {
					
					# Variable to switch to true if intersection is found:
					intersection_occurs <- FALSE
					
					# For each protected link:
					for(i in 1:nrow(protected_links)) {
						
						# Find if there are any intersections:
						intersections <- ArcIntersection(longitudes[match(as.character(protected_links[i, ])[1], colnames(intercontinent_links))], latitudes[match(as.character(protected_links[i, ])[1], colnames(intercontinent_links))], longitudes[match(as.character(protected_links[i, ])[2], colnames(intercontinent_links))], latitudes[match(as.character(protected_links[i, ])[2], colnames(intercontinent_links))], cut_point_1$long, cut_point_1$lat, cut_point_2$long, cut_point_2$lat, type = c("arc", "GC"), EarthRad = EarthRad)
						
						# If there are intersections update intersection_occurs:
						if(nrow(intersections) > 0) intersection_occurs <- TRUE
						
					}
					
					# Set up counter (to be used to warn user if loop gets stuck):
					counter <- 1
					
					# Whilst there is an intersection between the cut line and a protected link:
					while(intersection_occurs) {
						
						# Get new cut bearing (will be used to describe Great Circle made by cut):
						cut_bearing <- runif(1, 0, 360)
						
						# Pick new second point on Great Circle to use to describe it:
						cut_point_2 <- EndPoint(cut_point_1$longitude, cut_point_1$latitude, cut_bearing, min_separation, EarthRad = EarthRad)[c("longitude", "latitude")]
						
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
				cut_line_pole_1 <- EndPoint(cut_point_1$longitude, cut_point_1$latitude, (cut_bearing + 90) %% 360, 0.5 * pi * EarthRad, EarthRad = EarthRad)[c("longitude", "latitude")]
				
				# Get second pole to the cut line equator:
				cut_line_pole_2 <- EndPoint(cut_point_1$longitude, cut_point_1$latitude, (cut_bearing - 90) %% 360, 0.5 * pi * EarthRad, EarthRad = EarthRad)[c("longitude", "latitude")]
				
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
				
			}
			
		}
			
	# Case if no unprotected links:
	} else {
		
		# Make single clump from all continents:
		clumps <- paste(continent_numbers, collapse="&")
		
		# Warn user that no separation can be made:
		if(Warn) print("WARNING: Cannot separate continent as all links are protected.")
		
	}
	
	# Output clumps:
	return(clumps)
	
}
