#' Maps continents on Mercator projection
#'
#' Given a set of continent centres and their radius plots them on a Mercator projection
#' @param continent_centres A two-column matrix of decimalised longitudes (first column) and latitudes (second column).
#' @param radius A single radius value for each continent in kilometres.
#' @param xlim A single radius value for each continent in kilometres.
#' @param ylim A single radius value for each continent in kilometres.
#' @param resolution The number of points to use to plot each continent.
#' @param sea_colour The colour value to use to plot the sea.
#' @param land_colour The colour value used to plot the land.
#' @return A Mercator projection plot of each continent.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' MapContinents(continent_centres = matrix(c(-0, 0), ncol=2), radius = 2000)

MapContinents <- function(continent_centres, radius, xlim = c(-180, 180), ylim = c(-90, 90), resolution = 1000, sea_colour = "blue", land_colour = "green", ...) {
	
# Add grid overlay option:
# Add option to draw continent links to help visualise separation(s)/collision(s)
# Add ability to use other projection options than Mercator
# Ensure resolution has a minimum value of say, 100.
# Top-level conditionals to check values make sense, e.g., long-lat within limits.
# Allow borders of continents to be plotted?
	
	# Make empty plot:
	plot(x = c(1, 1), xlim = xlim, ylim = c(-90, 90), type = "n", xlab = "Longitude", ylab = "Latitude", ...)
	
	# Fill with blue background (sea):
	polygon(x = c(-180, 180, 180, -180), y = c(90, 90, -90, -90), border = NA, col = sea_colour)

	# Create empty list to store continent polygons:
	continent_polygons <- as.list(c(1:nrow(continent_centres)))
	
	# Create empty vector to store continents that cross the date line:
	continents_across_date_line <- vector(mode="numeric")
	
	# For each continent:
	for(i in 1:nrow(continent_centres)) {
	
		# Isolate continent centre:
		centre <- continent_centres[i, ]
		
		# Create empty matrix to store points that will describe the edge of the continent:
		edge_points <- matrix(nrow = 0, ncol = 2)
		
		# For as many points as are dictated by resolution get the jth point describing the edge of the continent:
		for(j in 1:resolution) edge_points <- rbind(edge_points, unlist(EndPoint(centre[1], centre[2], (360 / resolution) * j, radius)[c("longitude", "latitude")]))
		
		# Store full continent in list:
		continent_polygons[[i]] <- edge_points
		
		# Check to see if continent crosses date line and if so store the continent number:
		if(any(abs(diff(edge_points[, 1])) > 180)) continents_across_date_line <- c(continents_across_date_line, i)

		
	}

	# Get distances from North pole to each continents centre:
	North_pole_distances <- One2ManyGreatCircleDistance(0, 90, continent_centres[, 1], continent_centres[, 2])
	
	# Get distances from South pole to each continents centre:
	South_pole_distances <- One2ManyGreatCircleDistance(0, -90, continent_centres[, 1], continent_centres[, 2])

	# Get number(s) of continent(s) that sit over the North pole:
	continents_at_North_pole <- which((as.numeric(apply(matrix(North_pole_distances, nrow = 1), 2, all.equal, current = radius) != TRUE) + as.numeric(North_pole_distances < radius)) == 2)
	
	# Get number(s) of continent(s) that sit over the South pole:
	continents_at_South_pole <- which((as.numeric(apply(matrix(South_pole_distances, nrow = 1), 2, all.equal, current = radius) != TRUE) + as.numeric(South_pole_distances < radius)) == 2)
	
	# Remove continents that only cross the date line because they also overlap the pole:
	continents_across_date_line <- setdiff(continents_across_date_line, c(continents_at_North_pole, continents_at_South_pole))
	
	# If there are continent(s) at the North pole:
	if(length(continents_at_North_pole) > 0) {
	
		# For each continent at the North pole:
		for(i in continents_at_North_pole) {
			
			# Get signs (positive/negative) of longitudes (i.e., corresponding to East/West) of points in polygon:
			longitude_signs <- continent_polygons[[i]][, 1] / abs(continent_polygons[[i]][, 1])
			
			# Get breaks (where switch from positive to negative occurs):
			breaks <- which(abs(diff(longitude_signs)) > 0)
			
			# Reorder points using first break:
			continent_polygons[[i]] <- continent_polygons[[i]][c((breaks[1] + 1):resolution, 1:breaks[1]), ]
			
			# Build mini-matrix of additonal points (to take continent to date line and corners):
			extra_points <- matrix(c(c(-180, continent_polygons[[i]][which.min(continent_polygons[[i]][, 1]), 2]), c(-180, 90), c(180, 90), c(180, continent_polygons[[i]][which.max(continent_polygons[[i]][, 1]), 2])), ncol = 2, byrow = TRUE)
			
			# Add extra points in right part of continent to plot:
			continent_polygons[[i]] <- rbind(continent_polygons[[i]][1:which.max(continent_polygons[[i]][, 1]), ], extra_points, continent_polygons[[i]][setdiff(1:resolution, 1:which.max(continent_polygons[[i]][, 1])), ])
			
		}
	
	}
	
	# If there are continent(s) at the South pole:
	if(length(continents_at_South_pole) > 0) {

		# For each continent at the South pole:
		for(i in continents_at_South_pole) {

			# Get signs (positive/negative) of longitudes (i.e., corresponding to East/West) of points in polygon:
			longitude_signs <- continent_polygons[[i]][, 1] / abs(continent_polygons[[i]][, 1])
			
			# Get breaks (where switch from positive to negative occurs):
			breaks <- which(abs(diff(longitude_signs)) > 0)
			
			# Reorder points using first break:
			continent_polygons[[i]] <- continent_polygons[[i]][c((breaks[1] + 1):resolution, 1:breaks[1]), ]
			
			# Build mini-matrix of additonal points (to take continent to date line and corners):
			extra_points <- matrix(c(c(180, continent_polygons[[i]][which.max(continent_polygons[[i]][, 1]), 2]), c(180, -90), c(-180, -90), c(-180, continent_polygons[[i]][which.min(continent_polygons[[i]][, 1]), 2])), ncol = 2, byrow = TRUE)
			
			# Add extra points in right part of continent to plot:
			continent_polygons[[i]] <- rbind(continent_polygons[[i]][1:which.max(continent_polygons[[i]][, 1]), ], extra_points, continent_polygons[[i]][setdiff(1:resolution, 1:which.max(continent_polygons[[i]][, 1])), ])
	
		}
	
	}
	
	# If there are continent(s) at the North pole:
	if(length(continents_across_date_line) > 0) {
	
		# For each continent that crosses the date line:
		for(i in continents_across_date_line) {
		
			# Get breaks (where latitudinal "jump" occurs as crosses date line):
			breaks <- which(abs(diff(c(continent_polygons[[i]][, 1], continent_polygons[[i]][, 1]))) > 180)[1:2]
			
			# Get first side (first "half" of the polygon):
			side_one <- continent_polygons[[i]][(breaks[1] + 1):breaks[2], ]
			
			# Get second side (second "half" of the polygon):
			side_two <- continent_polygons[[i]][c(min(c(breaks[2] + 1, resolution)):resolution, 1:breaks[1]), ]
			
			# Add extra points to draw first polygon up to the date line:
			side_one <- rbind(side_one, matrix(c(as.vector(rep((side_one[1, 1] / abs(side_one[1, 1])) * 180, 2)), as.vector(c(side_one[nrow(side_one), 2], side_one[1, 2]))), ncol = 2))
			
			# Add extra points to draw second polygon up to the date line:
			side_two <- rbind(side_two, matrix(c(as.vector(rep((side_two[1, 1] / abs(side_two[1, 1])) * 180, 2)), as.vector(c(side_two[nrow(side_two), 2], side_two[1, 2]))), ncol = 2))
			
			# Add two sides as separate polygons to the end of the list:
			continent_polygons <- c(continent_polygons, list(side_one), list(side_two))
			
		}
		
		# Remove old polygon(s) from the list:
		continent_polygons <- continent_polygons[-continents_across_date_line]
	
	}
	
	# Plot continental polygons:
	for(i in 1:length(continent_polygons)) polygon(x = continent_polygons[[i]][, 1], y = continent_polygons[[i]][, 2], border = NA, col = land_colour)
	
}
