#' Global function to simulate evolution on a sphere
#'
#' Global function to generate simulated continents and clades on a sphere
#' @param N_steps The number of time steps in the simulation.
#' @param organism_multiplier The number of organism time steps to be taken per continent time step.
#' @param N_continents The (maximum) number of individual continents.
#' @param radius The radius of each circular continent.
#' @param start_configuration One of "random separate", "random overlap", "supercontinent", or "max separate".
#' @param squishiness A value from 0 (continents can never overlap) to 1 (continents can overlap completely).
#' @param stickiness Probability of two conjoined continents remaining conjoined in the next time step.
#' @param continent_speed_mean Mean continent speed (kilometers per time step).
#' @param continent_speed_sd Standard deviation of continent speed (kilometers per time step).
#' @param organism_step_sd Standard deviation used for random walk draws for organisms.
#' @param organism_step_shape Shape parameter to use for drawing bearings for organisms.
#' @param sweepstakes The probability of a sweepstakes dispersal occuring each time an organism hits the coast.
#' @param b Per-lineage birth (speciation) rate.
#' @param d Per-lineage death (extinction) rate.
#' @param EarthRad Radius of the Earth in kilometres.
#' @param polar Whether continents are allowed to exist exactly at the poles.
#' @return A magic table of awesomeness.
#' @details Nothing yet.
#' @author Laura C. Soul \email{lauracsoul@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @import ape
#' @import geiger
#' @examples
#' #DispersalSimulator(N_steps = 100, organism_multiplier = 5, N_continents = 2, radius = 2000,
#' #  start_configuration = "random separate", squishiness = 0.25, stickiness = 0.95,
#' #  continent_speed_mean = 5, continent_speed_sd = 2, organism_step_sd = 100, sweepstakes = 0, 
#' #  b = 0.1, d = 0.05, EarthRad = 6367.4447,  polar = FALSE)

# Total land area

# Use this line for debugging (sets values for input variables):
#N_steps = 200; organism_multiplier = 1; N_continents = 7; radius = 2000; start_configuration = "supercontinent"; squishiness = 0.25; stickiness = 0; continent_speed_mean = 5; continent_speed_sd = 2; organism_step_sd = 100; organism_step_shape = 1; b = 0.1; d = 0.05; sweepstakes = 0; EarthRad = 6367.4447; polar = FALSE

DispersalSimulator <- function(N_steps = 1000, organism_multiplier = 5, N_continents = 7, radius = 2000, start_configuration = "random separate", squishiness = 0.25, stickiness = 0.95, continent_speed_mean = 5, continent_speed_sd = 2, organism_step_sd = 100, organism_step_shape = 1, b = 0.1, d = 0.05, sweepstakes = 0, EarthRad = 6367.4447, polar = FALSE) {

# Need more top-level conditionals, e.g. N steps must be a positive integer, speed mean and sd must be non-negative
# Others may be caught by subfunctions so no need to repeat
# Add progress bar?
	
	# Start by picking continent start points:
	continent_starting_points <- StartingPoints(N_continents = N_continents, radius = radius, start_configuration = start_configuration, squishiness = squishiness, EarthRad = EarthRad, polar = polar)
	
# Are any continents joined at start?:

	# Get minimum_continental separation:
	min_separation <- (1 - squishiness) * radius * 2

	# Get list of separate continents:
	separate_continents <- HowManySeparateContinents(min_separation, continent_starting_points[, "Longitude"], continent_starting_points[, "Latitude"])

	# Get list of touching continents (to be used later for whether dispersal is allowable or not):
	touching_continents <- HowManySeparateContinents((radius * 2), continent_starting_points[, "Longitude"], continent_starting_points[, "Latitude"])
	
# Assign Euler poles to each separate continent:
# THESE MAY NOT BE EQUALLY DISTRIBUTED ACROSS THE PLANET, MAYBE RECONSIDER HOW TO DRAW THESE
	
	# Randomly draw longitudes for each separated continent:
	euler_pole_longitudes <- runif(length(separate_continents), -180, 180)
	
	# Randomly draw latitudes for each separated continent:
	euler_pole_latitudes <- runif(length(separate_continents), -90, 90)
	
	# Ensure no latitude is directly at a pole (North or South) by redrawing if any are found (causes an issue with bearings - e.g., what is 0-degrees from North pole):
	if((sum(euler_pole_latitudes == 90) + sum(euler_pole_latitudes == -90)) > 0) euler_pole_latitudes <- runif(length(separate_continents), -90, 90)

# Assign speeds to each separate continent:
	
	# Create empty vector to store degrees per step (effectively the speed of movement) for each continent:
	degrees_per_step <- vector(mode="numeric")
	
	# For each separate continent:
	for(i in 1:length(separate_continents)) {
		
		# Get Greate Circle distances from Euler pole to each continent centre:
		euler_GC_distances <- One2ManyGreatCircleDistance(euler_pole_longitudes[i], euler_pole_latitudes[i], continent_starting_points[as.numeric(unlist(strsplit(separate_continents[i], "&"))), "Longitude"], continent_starting_points[as.numeric(unlist(strsplit(separate_continents[i], "&"))), "Latitude"])
		
		# Find GC distance to furthest continent (closest to euler pole "equator") in cluster (as speed will be assigned based on this):
		furthest_continent_GC_distance <- (0.5 * pi * EarthRad) - min(abs(euler_GC_distances - (0.5 * pi * EarthRad)))
		
		# Randomly draw a continent speed:
		continent_speed <- rnorm(n = 1, mean = continent_speed_mean, sd = continent_speed_sd)
		
		# If a negative speed is picked then redraw:
		while(continent_speed < 0) continent_speed <- rnorm(n = 1, mean = continent_speed_mean, sd = continent_speed_sd)
		
		# Set degree change per step (effectively the speed):
		degrees_per_step[i] <- continent_speed / (2 * pi * furthest_continent_GC_distance) * 360
		
	}

# Starting information for the tree:
	
# MAYBE DISTINGUISH HERE BETWEEN REGULAR AND SWEEPSTAKES DISPERSALS	
	# Create marix that will store each intercontinental dispersal:
	dispersals <- matrix(nrow = 0, ncol = 5, dimnames = list(c(), c("From_continent", "To_continent", "Branch", "Continent_time_step", "Animal_time_step")))
	
	# Pick a random continent on which the clade begins:
	begin_cont <- ceiling(runif(1, 0, N_continents))
	
	# Start by drawing two values from long-lat limits as though continent were centered at equator-Greenwich meridean intersection (means degrees are same distance value):
	begin_coordinates <- runif(n = 2, min = -radius / (2 * pi * EarthRad) * 360, max = radius / (2 * pi * EarthRad) * 360)
	
	# Whilst the draw is not on the continent redraw long and lat values:
	while(GreatCircleDistanceFromLongLat(0, 0, begin_coordinates[1], begin_coordinates[2]) > radius) begin_coordinates <- runif(n = 2, min = -radius / (2 * pi * EarthRad) * 360, max = radius / (2 * pi * EarthRad) * 360)
	
	# Get bearing from continent centre to point where clade begins:
	bearing_to_clade_start <- BearingBetweenTwoLongLatPoints(0, 0, begin_coordinates[1], begin_coordinates[2])

	# Get distance from continent centre to point where clade begins:
	distance_to_clade_start <- GreatCircleDistanceFromLongLat(0, 0, begin_coordinates[1], begin_coordinates[2], EarthRad = EarthRad, Warn = FALSE)
	
	# Get coordinates of point where clade begins:
	life_begins <- EndPoint(continent_starting_points[begin_cont, "Longitude"], continent_starting_points[begin_cont, "Latitude"], bearing_to_clade_start, distance_to_clade_start, EarthRad = EarthRad)
	
	extra.rows <- matrix(NA, nrow = 2, ncol = (N_steps * organism_multiplier) + 1)
	organism_lat_matrix <- matrix(nrow = 2, ncol = (N_steps * organism_multiplier) + 1)
    organism_long_matrix <- matrix(nrow = 2, ncol = (N_steps * organism_multiplier) + 1)
    organism_lat_matrix[, 1] <- life_begins$latitude
    organism_long_matrix[, 1] <- life_begins$longitude
    rownames(organism_lat_matrix) <- rep(begin_cont, 2)
    rownames(organism_long_matrix) <- rep(begin_cont, 2)
    edge <- rbind(c(1, 2), c(1, 3)) # this is a starting edge matrix
    edge.length <- rep(NA, 2)
    stem.depth <- numeric(2)
    alive <- rep(TRUE, 2) # marker for live lineages
    ot <- 0 # time(step) at any point in the tree
    next.node <- 4

# Moving continents:
	
	# Create empty lists to store which circles are in each supercontinent or that are touching allowing dispersal (new elements will be added only when these change):
	touching <- linked <- list()
	
	# Use starting values to fill first item in list:
	linked[[1]] <- separate_continents

	# Use starting values to fill first item in list:
	touching[[1]] <- touching_continents
	
	# Number based on first time step (technically this should be zero as no steps have occurred, but here we are using one):
	names(linked)[[1]] <- "1:1"
	
	# Number based on first time step (technically this should be zero as no steps have occurred, but here we are using one):
	names(touching)[[1]] <- "1:1"
	
	# List to store which circles are in each supercontinent (new element added only when it changes)

	# Array to store the positions of every continent at each time step
	position <- array(NA, c(N_continents, N_steps + 1, 2), c("continent", "timestep", "coordinate"))
	position[, 1, 1] <- continent_starting_points[, "Longitude"]
	position[, 1, 2] <- continent_starting_points[, "Latitude"]

	temp_position <- matrix(nrow = N_continents, ncol = 2)

# REPLACE THIS WITH PROGRESS BAR AT SOME POINT
	progress <- seq(1, N_steps, floor(N_steps / 50))
	cat("Starting time steps \n")
	cat("Progress \n")
	cat("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")

	# Main continent time step loop:
	for (t in 2:(N_steps + 1)) {
		
		# Distances continents are apart before moving:
		starting_distances <- GreatCircleDistanceMatrix(position[, t - 1, 1], position[, t - 1, 2])
		
		for (k in 1:N_continents) {
			
			# Find current longlat of the kth continent:
			start_long <- position[k, t - 1, 1]
			start_lat <- position[k, t - 1, 2]

			# Identify which supercontinent, and therefore which element of the euler pole and speed vectors, the circle k belongs to
			where <- WhichSupercontinent(k, tail(linked, n = 1)[[1]])

			# Find distance of circle from pole
			distance <- GreatCircleDistanceFromLongLat(long1 = start_long,lat1=start_lat, long2 = euler_pole_longitudes[where], lat2 = euler_pole_latitudes[where], Warn = FALSE)
			
			# Find bearing of circle from pole
			init_bearing <- BearingBetweenTwoLongLatPoints(euler_pole_longitudes[where], euler_pole_latitudes[where], start_long, start_lat)

			# Find the new bearing of circle from pole according to the speed specified
			new_bearing <- (init_bearing + degrees_per_step[where]) %% 360

			# Find the new location of the circle
			new_loc <- EndPoint(euler_pole_longitudes[where], euler_pole_latitudes[where], new_bearing, distance, EarthRad = EarthRad)

			# Add the new loction to the position matrix
			temp_position[k, 1] <- new_loc$longitude
			temp_position[k, 2] <- new_loc$latitude

		}

		# Matrix of the distances between each continent in their new positions, before these are set
		new_distances <- GreatCircleDistanceMatrix(temp_position[, 1], temp_position[, 2])
		
		# Matrix of euler poles and speeds for each individual continent for moving the animals later
		organism_mover <- matrix(nrow = N_continents, ncol = 3)
		for (clump in 1:length(tail(linked, n = 1)[[1]])) {
			which_conts <- as.numeric(strsplit(tail(linked, n = 1)[[1]][[clump]], "&")[[1]])
			organism_mover[which_conts, 1] <- rep(euler_pole_longitudes[clump], length(which_conts))
			organism_mover[which_conts, 2] <- rep(euler_pole_latitudes[clump], length(which_conts))
			# This bit will get overwritten if there are collisions
			organism_mover[which_conts, 3] <- rep(degrees_per_step[clump], length(which_conts))
		}
		
		# Making a vector of whether any continents have collided and thus should potentially be linked:
		collisions <- matrix(nrow = 0, ncol = 2)
		for (q in 1:(N_continents - 1)) {
			for (p in (q + 1):N_continents) {
				comp1 <- WhichSupercontinent(p, tail(linked, n = 1)[[1]]) == WhichSupercontinent(q, tail(linked, n = 1)[[1]])
				comp2 <- all.equal(new_distances[p, q], min_separation) == TRUE || new_distances[p, q] < min_separation
				if (comp1 != TRUE && comp2 == TRUE) {
					where_1 <- WhichSupercontinent(p, tail(linked, n = 1)[[1]])
					where_2 <- WhichSupercontinent(q, tail(linked, n = 1)[[1]])
					continent_1_euler_longitude <- organism_mover[p, 1]
					continent_1_euler_latitude <- organism_mover[p, 2]
					continent_2_euler_longitude <- organism_mover[q, 1]
					continent_2_euler_latitude <- organism_mover[q, 2]
					continent_1_degrees_per_step <- organism_mover[p, 3]
					continent_2_degrees_per_step <- organism_mover[q, 3]
					proportion_to_reverse <- ColliderReverser(min_separation, position[p, t - 1, 1], position[p, t - 1, 2], temp_position[p, 1], temp_position[p, 2], position[q, t - 1, 1], position[q, t - 1, 2], temp_position[q, 1], temp_position[q, 2], continent_1_euler_longitude, continent_1_euler_latitude, continent_2_euler_longitude, continent_2_euler_latitude, continent_1_degrees_per_step, continent_2_degrees_per_step, EarthRad = EarthRad, Warn = FALSE)
					collisions <- rbind(collisions, sort(c(p, q)))
					rownames(collisions)[nrow(collisions)] <- paste(c(sort(c(tail(linked, n = 1)[[1]][where_1], tail(linked, n = 1)[[1]][where_2])), proportion_to_reverse), collapse = " ")
				}
			}
		}

		# Matrices to store permanent collisions:
		perm_collisions <- clump_collisions <- matrix(nrow = 0, ncol = 2)
		
		# Move continents until all collisions have been accounted for:
		while(nrow(collisions) > 0) {
			
			proportions <- as.numeric(unlist(lapply(strsplit(rownames(collisions), " "), '[', 3)))
			
			# Which continents collide first?:
			first_collisions <- which(proportions == min(proportions))

			# What other collisions occur?:
			other_collisions <- setdiff(1:nrow(collisions), first_collisions)
			
			# Find the continents involved in the earliest collision(s):
			cont_involved <- as.numeric(unique(unlist(strsplit(unlist(lapply(strsplit(rownames(collisions)[first_collisions], " "), '[', 1:2)), "&"))))

			# If there are other collisions (need to check if they include same continents):
			if(length(other_collisions) > 0) {
				
				# Get list of continents involved in each other (later) collision(s):
				other_cont_involved <- lapply(lapply(lapply(lapply(strsplit(rownames(collisions)[other_collisions], " "), '[', 1:2), strsplit, split = "&"), unlist), as.numeric)
			
				# Add any continents attached to those involved in first collision(s) (as these should be moved back too):
				cont_involved <- unique(c(cont_involved, unlist(other_cont_involved[which(unlist(lapply(lapply(lapply(other_cont_involved, match, cont_involved), sort), length)) > 0)])))

			}
				
			# Update temporary positions based on new change in bearing for all the continents the other clump
			for (rev in cont_involved) {
				
# Make this faster by doing by clump instead of by continent?
				
				start_long <- position[rev, t - 1, 1]
				start_lat <- position[rev, t - 1, 2]
				min_proportion <- min(proportions)
				if(min_proportion != 0) {
					new_bearing <- init_bearing <- BearingBetweenTwoLongLatPoints(organism_mover[rev, 1], organism_mover[rev, 2], start_long, start_lat)
					distance <- GreatCircleDistanceFromLongLat(long1 = start_long, lat1 = start_lat, long2 = organism_mover[rev, 1], lat2 = organism_mover[rev, 2], Warn = FALSE)
					
					# Store modified degrees per step:
					organism_mover[rev, 3] <- (min_proportion * organism_mover[rev, 3])
					new_bearing <- (init_bearing + organism_mover[rev, 3]) %% 360
					new_loc <- EndPoint(organism_mover[rev, 1], organism_mover[rev, 2], new_bearing, distance, EarthRad = EarthRad)
					temp_position[rev, 1] <- new_loc$longitude
					temp_position[rev, 2] <- new_loc$latitude
				} else {
					temp_position[rev, 1] <- start_long
					temp_position[rev, 2] <- start_lat
					organism_mover[rev, 3] <- 0
				}
			}

			# Add to matrix of definite collisions that cannot be separated in the next step:
			clump_collisions <- rbind(clump_collisions, matrix(unlist(lapply(strsplit(rownames(collisions)[first_collisions], " "), '[', 1:2)), ncol = 2, byrow = TRUE))
			
			
			perm_collisions <- rbind(perm_collisions, collisions[first_collisions, ])
			
			# Recalculate the new distances:
			new_distances <- GreatCircleDistanceMatrix(temp_position[, 1], temp_position[, 2])
		
			# Check whether any other collisions have still occurred
			collisions <- matrix(nrow = 0, ncol = 2)
			
			for (q in 1:(N_continents - 1)) {
				
				for (p in (q + 1):N_continents) {
					
					comp1 <- WhichSupercontinent(p, tail(linked, n = 1)[[1]]) == WhichSupercontinent(q, tail(linked, n = 1)[[1]])
					comp2 <- all.equal(new_distances[p, q], min_separation) == TRUE || new_distances[p, q] < min_separation
					
					if (comp1 != TRUE && comp2 == TRUE) {
						where_1 <- tail(linked, n = 1)[[1]][WhichSupercontinent(p, tail(linked, n = 1)[[1]])]
						where_2 <- tail(linked, n = 1)[[1]][WhichSupercontinent(q, tail(linked, n = 1)[[1]])]
						
						# If collision has not already been recorded:
						if(!any(apply(matrix(as.numeric(clump_collisions == where_1) + as.numeric(clump_collisions == where_2), ncol = 2), 1, sum) == 2)) {
						
							where_1 <- WhichSupercontinent(p, tail(linked, n = 1)[[1]])
							where_2 <- WhichSupercontinent(q, tail(linked, n = 1)[[1]])
							continent_1_euler_longitude <- organism_mover[p, 1]
							continent_1_euler_latitude <- organism_mover[p, 2]
							continent_2_euler_longitude <- organism_mover[q, 1]
							continent_2_euler_latitude <- organism_mover[q, 2]
							continent_1_degrees_per_step <- organism_mover[p, 3]
							continent_2_degrees_per_step <- organism_mover[q, 3]
							proportion_to_reverse <- ColliderReverser(min_separation, position[p, t - 1, 1], position[p, t - 1, 2], temp_position[p, 1], temp_position[p, 2], position[q, t - 1, 1], position[q, t - 1, 2], temp_position[q, 1], temp_position[q, 2], continent_1_euler_longitude, continent_1_euler_latitude, continent_2_euler_longitude, continent_2_euler_latitude, continent_1_degrees_per_step, continent_2_degrees_per_step, EarthRad = EarthRad, Warn = FALSE)
							collisions <- rbind(collisions, sort(c(p, q)))
							rownames(collisions)[nrow(collisions)] <- paste(c(sort(c(tail(linked, n = 1)[[1]][where_1], tail(linked, n = 1)[[1]][where_2])), proportion_to_reverse), collapse = " ")
							
						}
						
					}
					
				}
				
			}
			
		}

		# When it finishes while loop can now 'ossify' the positions and move to selecting separations
		position[, t, 1] <- temp_position[, 1]
		position[, t, 2] <- temp_position[, 2]

# Finding out if anything gets separated
		
		# How many separate continents are there now? (after collisions may have occurred):
		separate_continents <- HowManySeparateContinents(min_separation, position[, t, 1], position[, t, 2])
		
		# Make random uniform draws for each continental cluster:
		separation_draws <- runif(length(grep("&", separate_continents)))
		
		# Case if a separation occurs (drawn value is equal to or exceeds stickiness):
		if(any(separation_draws >= stickiness)) {
		
			# Get clusters of continents to split apart:
			clusters_to_split <- separate_continents[grep("&", separate_continents)][which(separation_draws >= stickiness)]
			
			# For each cluster to split apart:
			for(i in 1:length(clusters_to_split)) {
				
				# Get numbers of continents involved:
				cluster_continent_numbers <- sort(strsplit(clusters_to_split[i], "&")[[1]])
				
				# Get longitudes of continents in cluster:
				cluster_longitudes <- position[as.numeric(cluster_continent_numbers), t, 1]

				# Get latitudes of continents in cluster:
				cluster_latitudes <- position[as.numeric(cluster_continent_numbers), t, 2]
				
				# Make protected links an empty matrix:
				cluster_protected_links <- matrix(nrow = 0, ncol = 2)
				
				# If there are recent collisions (that will potentially need to be excluded from new splits):
				if(length(perm_collisions) > 0) {
					
					# For each new collision:
					for(j in 1:nrow(perm_collisions)) {
						
						# If new collision occurs within the cluster:
						if(length(intersect(as.character(perm_collisions[j, ]), cluster_continent_numbers)) == 2) {
						
							# Add new collision to protected links list:
							cluster_protected_links <- rbind(cluster_protected_links, perm_collisions[j, ])
							
						}
						
					}

				}
				
				# Get splits (new clusters) of separated cluster:
				splits <- ContinentSplitter(min_separation, cluster_longitudes, cluster_latitudes, cluster_continent_numbers, cluster_protected_links, EarthRad, Warn = FALSE)
				
				# Update separate continents vector:
				separate_continents <- sort(c(separate_continents[-match(clusters_to_split[i], separate_continents)], splits))
				
			}
			
		}
		
		# Get current list of touching continents (to be used later for whether dispersal is allowable or not):
		touching_continents <- HowManySeparateContinents((radius * 2), position[, t, 1], position[, t, 2])
		
		# If touching relationships between continetns have changed:
		if (paste(sort(touching_continents), collapse="") != paste(sort(tail(touching, n = 1)[[1]]), collapse="")) {
			
			# Add new element to the touching list:
			touching <- c(touching, list(touching_continents))

			# Update time step for previous continental configuration:
			names(touching)[(length(touching) - 1)] <- paste(strsplit(names(touching[(length(touching) - 1)]), ":")[[1]][1], ":", t - 1, sep="")
			
			# Add time step to new continental configuration:
			names(touching)[length(touching)] <- paste(t, ":", t, sep="")

		}

# Now change the Euler poles and speeds

		# If the continental clustering has changed:
		if (paste(sort(separate_continents), collapse="") != paste(sort(tail(linked, n = 1)[[1]]), collapse="")) {
			
			# Add new continental configuration to linked list
			linked <- c(linked, list(separate_continents))

			# Update time step for previous continental configuration:
			names(linked)[(length(linked) - 1)] <- paste(strsplit(names(linked[(length(linked) - 1)]), ":")[[1]][1], ":", t - 1, sep="")
			
			# Add time step to new continental configuration:
			names(linked)[length(linked)] <- paste(t, ":", t, sep="")
			
			# Select continents that are different to previous time step
			toKeep <- c(tail(linked, n = 1)[[1]], tail(linked, n = 2)[[1]])[duplicated(c(tail(linked, n = 1)[[1]], tail(linked, n = 2)[[1]]))]

			# Vectors for new poles and speeds
			new_euler_latitudes <- rep(NA, length(separate_continents))
			new_euler_longitudes <- rep(NA, length(separate_continents))
			new_degrees_per_step <- rep(NA, length(separate_continents))

			if (length(toKeep) != 0){
				# Inherit previous poles and speeds for those clumps that remain the same
				for (y in 1:length(toKeep)) {
					new_euler_longitudes[match(toKeep[y], separate_continents)] <- euler_pole_longitudes[match(toKeep[y], tail(linked, n = 2)[[1]])]
					new_euler_latitudes[match(toKeep[y], separate_continents)] <- euler_pole_latitudes[match(toKeep[y], tail(linked, n = 2)[[1]])]
					new_degrees_per_step[match(toKeep[y], separate_continents)] <- degrees_per_step[match(toKeep[y], tail(linked, n = 2)[[1]])]
				}
			}

			# Make new Euler poles
			new_euler_longitudes[is.na(new_euler_longitudes)] <- runif((length(separate_continents) - length(toKeep)), -180, 180)

			changed_euler_latitudes <- runif((length(separate_continents) - length(toKeep)), -90, 90)
			while((sum(changed_euler_latitudes == 90) + sum(changed_euler_latitudes == -90)) > 0) changed_euler_latitudes <- runif((length(separate_continents) - length(toKeep)), -90, 90)
			
			new_euler_latitudes[is.na(new_euler_latitudes)] <- changed_euler_latitudes

			# Get Great Circle distances from Euler pole to each continent centre:
			for (l in 1:(length(separate_continents) - length(toKeep))) {
				
				# Find the first NA to fill in:
				changer <- match(NA, new_degrees_per_step)
				
				# Find the continents that are in the clump whose speed is going to change:
				rows <- as.numeric(unlist(strsplit(separate_continents[changer], "&")))
				
				# Find the new euler pole for that clump:
				euler_long <- new_euler_longitudes[changer]
				euler_lat <- new_euler_latitudes[changer]
				
				# Find the distances of each of the continents in the clump from their new euler pole:
				euler_GC_distances <- One2ManyGreatCircleDistance(euler_long, euler_lat, position[rows, t, 1], position[rows, t, 2])
				
				# Find which of those is closest to the equator:
				furthest_continent_GC_distance <- (0.5 * pi * EarthRad) - min(abs(euler_GC_distances - (0.5 * pi * EarthRad)))
	
				# Randomly draw a continent speed:
				continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
				# If a negative speed is picked then redraw:
				while(continent_speed < 0) continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
				# Set degree change per step (effectively the speed):
				new_degrees_per_step[changer] <- continent_speed / (2 * pi * furthest_continent_GC_distance) * 360
				
			}
		
			euler_pole_longitudes <- new_euler_longitudes
			euler_pole_latitudes <- new_euler_latitudes
			degrees_per_step <- new_degrees_per_step
		}

#Need to add in the tree bit
		
		# Move the animals with the continents
		for (cont in 1:N_continents) {
			
			# Find the rows that are in the continent of focus
			moving <- which(rownames(organism_long_matrix) == cont)

			if (length(moving) == 0) {
				next 
			} else {
				# Find out how many degrees around the euler pole that continent went
				organism_degrees <- organism_mover[cont, 3]
				# loops through all the organisms that are currently living
				for (organism in 1:length(moving)) {
					organism_row <- moving[organism]
					if (is.na(organism_long_matrix[organism_row, ot + 1])) {
						next
					} else {
						first_long <- organism_long_matrix[organism_row, ot + 1]
						first_lat <- organism_lat_matrix[organism_row, ot + 1]
						organism_distance <- GreatCircleDistanceFromLongLat(organism_mover[cont, 1], organism_mover[cont, 2], first_long, first_lat)
						organism_bearing <- BearingBetweenTwoLongLatPoints(organism_mover[cont, 1], organism_mover[cont, 2], first_long, first_lat)
						new_organism_bearing <- (organism_bearing + organism_degrees) %% 360
						new_organism_loc <- EndPoint(organism_mover[cont, 1], organism_mover[cont, 2], new_organism_bearing, organism_distance, EarthRad = EarthRad)
						organism_long_matrix[organism_row, ot + 1] <- new_organism_loc$longitude
    					organism_lat_matrix[organism_row, ot + 1] <- new_organism_loc$latitude
						
					}
					
				}
				
			}
			
		}

		for (f in 1:organism_multiplier) {
# PERHAPS ALLOW LOOP TO COMPLETE ANYWAY ONCE COMPLETE EXTINCTION OCCURS
            if (sum(alive) == 0) stop("METEOR IMPACT!! EVERYTHING HAS GONE EXTINCT! AAAAHHHHHH!!!")
            dt <- 1
            ot <- ot + dt
            for (m in 1:nrow(organism_lat_matrix)) {
                if (alive[m]) {
                    starting <- c(organism_long_matrix[m, ot], organism_lat_matrix[m, ot])
                    moveto <- EndPoint(starting[1], starting[2], RandomBearingDraw(n = 1, shape_parameter = organism_step_shape), abs(rnorm(1, 0, organism_step_sd)), EarthRad = EarthRad) #generates a random walk step and calculates new position
					on_cont <- as.numeric(rownames(organism_lat_matrix)[m])
                    friends <- as.numeric(strsplit(tail(touching, n = 1)[[1]][WhichSupercontinent(on_cont, tail(touching, n = 1)[[1]])], "&")[[1]])
                    
					# Get distance new move-to step will take organism away from the current continents centre:
					dist_from_centre <- GreatCircleDistanceFromLongLat(position[on_cont, t, 1], position[on_cont, t, 2], moveto$longitude, moveto$latitude)
                    
					# If other continents are in contact with the currently inhabited continent:
					if (length(friends) > 1) {
						
                    	all_dist <- vector(length = length(friends))
                    	for (g in 1:length(friends)) all_dist[g] <- GreatCircleDistanceFromLongLat(position[friends[g], t, 1], position[friends[g], t, 2], moveto$longitude, moveto$latitude)
						new_cont <- friends[which(all_dist == min(all_dist))[1]]
						
						# If the current move-to step takes the organism out to sea:
						if(all.equal(min(all_dist), radius) != TRUE && min(all_dist) > radius) {
							
							# Case if a sweeptakes dispersal occurs:
							if(runif(n = 1, min = 0, max = 1) <= sweepstakes) {
							
								# Find landing spot after sweepstakes dispersal:
								landing_spot <- MagicPortal(start_longitude = as.vector(starting[1]), start_latitude = as.vector(starting[2]), end_longitude = moveto$longitude, end_latitude = moveto$latitude, start_continent = on_cont, touching_continents = as.vector(unlist(tail(touching, n = 1))), continent_centres = position[,t,], continent_radius = radius, EarthRad = EarthRad)
								
								# Get new continent after sweepstakes dispersal:
								new_cont <- landing_spot$continent
								
								# Update the move-to variable with the landing spot for the sweepstakes dispersal:
								moveto <- landing_spot[c("longitude", "latitude")]
								
							# Case if no sweepstakes dispersal occurs:
							} else  {
								
								# Whilst the current move-to step takes the organism off the continent:
								while (all.equal(min(all_dist), radius) != TRUE && min(all_dist) > radius) {
									
									# Draw a new move-to step:
									moveto <- EndPoint(starting[1], starting[2], RandomBearingDraw(n = 1, shape_parameter = organism_step_shape), abs(rnorm(1, 0, organism_step_sd)), EarthRad = EarthRad)
									temp_all_dist <- vector(length = length(friends))
									for (g in 1:length(friends)) {
										temp_all_dist[g] <- GreatCircleDistanceFromLongLat(position[friends[g], t, 1], position[friends[g], t, 2], moveto$longitude, moveto$latitude)
									}
									all_dist <- temp_all_dist
									new_cont <- friends[which(all_dist == min(all_dist))][1]
									
								}
								
							}
							
						}
                    	 
                    	organism_long_matrix[m, ot + 1] <- moveto$longitude
                    	organism_lat_matrix[m, ot + 1] <- moveto$latitude
                    	rownames(organism_long_matrix)[m] <- new_cont
                    	rownames(organism_lat_matrix)[m] <- new_cont
                    	if(new_cont != on_cont) dispersals <- rbind(dispersals, c(on_cont, new_cont, m, t, ot + 1))
						
					# If currently inhabited continent is isolated:
                    } else {
						
						# If the current move-to step takes the organism out to sea:
						if (all.equal(as.vector(dist_from_centre), radius) != TRUE && dist_from_centre > radius) {
							
							# Case if a sweeptakes dispersal occurs:
							if(runif(n = 1, min = 0, max = 1) <= sweepstakes) {
								
								# Get landing spot following sweepstakes dispersal:
								landing_spot <- MagicPortal(start_longitude = as.vector(starting[1]), start_latitude = as.vector(starting[2]), end_longitude = moveto$longitude, end_latitude = moveto$latitude, start_continent = on_cont, touching_continents = as.vector(unlist(tail(touching, n = 1))), continent_centres = position[,t,], continent_radius = radius, EarthRad = EarthRad)
								
								# Update continent taxon found on:
								new_cont <- landing_spot$continent
								
								# Update move-to variable with new position:
								moveto <- landing_spot[c("longitude", "latitude")]
								
								# Update continent organism found on:
								rownames(organism_long_matrix)[m] <- new_cont
								
								# Update continent organism found on:
								rownames(organism_lat_matrix)[m] <- new_cont
								
								# If continent is different to starting continent:
								if(new_cont != on_cont) {
									
									# Add sweepstakes dispersal to dispersals matrix:
									dispersals <- rbind(dispersals, c(on_cont, new_cont, m, t, ot + 1))
									
								}
								
							# Case if no sweepstakes dispersal occurs:
							} else  {
								
								# Whilst the current move-to step take the organism off the continent:
								while (all.equal(as.vector(dist_from_centre), radius) != TRUE && dist_from_centre > radius) {
									
									# Draw a new move-to step:
									moveto <- EndPoint(starting[1], starting[2], RandomBearingDraw(n = 1, shape_parameter = organism_step_shape), abs(rnorm(1, 0, organism_step_sd)), EarthRad = EarthRad)
									
									# Update distance from center:
									dist_from_centre <- GreatCircleDistanceFromLongLat(position[on_cont, t, 1], position[on_cont, t, 2], moveto$longitude, moveto$latitude)
								
								}
								
							}
							
						}
						
						# Update organism position in longitude matrix:
                    	organism_long_matrix[m, ot + 1] <- moveto$longitude

						# Update organism position in latitude matrix:
                    	organism_lat_matrix[m, ot + 1] <- moveto$latitude
						
                	}
                }
            }

            r <- runif(1)
            if (r <= b / (b + d)) { ###4 #this creates a bifucation in the tree
                random_lineage <- round(runif(1, min = 1, max = sum(alive)))
                e <- matrix(edge[alive, ], ncol = 2)
                parent <- e[random_lineage, 2]
                x <- which(edge[, 2] == parent)
                new.rows.lat <- new.rows.long <- extra.rows
                new.rows.lat[, ot + 1] <- organism_lat_matrix[x, ot + 1] #updates the long and lat matricies with new coordiantes for each lineage
                new.rows.long[, ot + 1] <- organism_long_matrix[x, ot + 1]
                rownames(new.rows.lat) <- rep(rownames(organism_lat_matrix)[x], 2)
                rownames(new.rows.long) <- rep(rownames(organism_long_matrix)[x], 2)
                organism_lat_matrix <- rbind(organism_lat_matrix, new.rows.lat)
                organism_long_matrix <- rbind(organism_long_matrix, new.rows.long)
                alive[alive][random_lineage] <- FALSE
                edge <- rbind(edge, c(parent, next.node), c(parent, next.node + 1))
                next.node <- next.node + 2
                alive <- c(alive, TRUE, TRUE)
                stem.depth <- c(stem.depth, ot, ot)
                edge.length[x] <- ot - stem.depth[x]
                edge.length <- c(edge.length, NA, NA)
            } else {###4 This terminates one of the current lineages on the tree
                random_lineage <- round(runif(1, min = 1, max = sum(alive)))
                edge.length[alive][random_lineage] <- ot - stem.depth[alive][random_lineage]
                alive[alive][random_lineage] <- FALSE
            }###4
        }
        if (any(progress == t)) {
        	cat("- ")
        }
    }

    # Things to do right at the end to turn it into a phylo object
    edge.length[alive] <- ot - stem.depth[alive]
    n <- -1
    for (f in 1:max(edge)) {
        if (any(edge[, 1] == f)) {
           	edge[which(edge[, 1] == f), 1] <- n
           	edge[which(edge[, 2] == f), 2] <- n
           	n <- n - 1
        }
    }
    
    edge[edge > 0] <- 1:sum(edge > 0)
   	tip.label <- 1:sum(edge > 0)
    mode(edge) <- "character"
    mode(tip.label) <- "character"
    obj <- list(edge = edge, edge.length = edge.length, tip.label=tip.label)
    class(obj) <- "phylo"
    obj <- old2new.phylo(obj)
    obj$tip.label = paste("s", 1:Ntip(obj), sep = "")
    
	# Add final time step to continental configurations:
	names(linked)[length(linked)] <- paste(strsplit(names(linked)[length(linked)], ":")[[1]][1], ":", t, sep="")
	
	# Add final time step to continental configurations:
	names(touching)[length(touching)] <- paste(strsplit(names(touching)[length(touching)], ":")[[1]][1], ":", t, sep="")

	output <- list(position, linked, touching, obj, organism_long_matrix, organism_lat_matrix, dispersals)
	
	names(output) <- c("continent_positions", "continent_clusters", "continent_overlaps", "tree", "organism_longitudes", "organism_latitudes", "dispersals")
	
	return(output)
	
}
