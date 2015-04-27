#' Global function to simulate evolution on a sphere
#'
#' Global function to generate simulated continents and clades on a sphere
#'
#' @param N_steps The number of time steps in the simulation.
#' @param N_continents The (maximum) number of individual continents.
#' @param radius The radius of each circular continent.
#' @param start_configuration One of "random separate", "random overlap", "supercontinent", or "max separate".
#' @param squishiness A value from 0 (continents can never overlap) to 1 (continents can overlap completely).
#' @param stickiness Probability of two conjoined continents remaining conjoined in the next time step.
#' @param continent_speed_mean Mean continent speed (kilometers per time step).
#' @param continent_speed_sd Standard deviation of continent speed (kilometers per time step).
#' @param EarthRad Radius of the Earth in kilometres.
#' @param polar Whether continents are allowed to exist exactly at the poles.
#' @return A magic table of awesomeness.
#' @details Nothing yet.
#'
#' @examples
#' EverythingFunction(N_steps = 10, N_continents = 2, radius = 2000,
#'   start_configuration = "supercontinent", squishiness = 0.25, stickiness = 0.95,
#'   continent_speed_mean = 5, continent_speed_sd = 2, EarthRad = 6367.4447,
#'   polar = FALSE)

# Outputs:
# - N separated continents (distance matrix with values less than 2 radii)
# - Long-lat of each circle centre
# - Bearings after each step change
# - Total land area all circles - minus

EverythingFunction <- function(N_steps = 1000, organism_multiplier = 1, N_continents = 7, radius = 2000, start_configuration = "supercontinent", squishiness = 0.25, stickiness = 0.95, continent_speed_mean = 5, continent_speed_sd = 2, EarthRad = 6367.4447, polar = FALSE) {

# Need more top-level conditionals, e.g. N steps must be a positive integer, speed mean and sd must also be positive
# Others may be caught by subfunctions so no need to repeat
# Add progress bar?
	
	# Subfunction to find out which supercontinent a circle belongs to:
	which_supercontinent <- function(cont, sprconts){
		cont <- as.character(cont)
		result <- which(unlist(lapply(lapply(lapply(strsplit(sprconts, "&"), match, cont), sort), length)) == 1)
		return(result)
	}
	
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
	
	# Randomly draw longitudes for each separated continent:
	euler_pole_longitudes <- runif(length(separate_continents), -180, 180)
	
	# Randomly draw latitudes for each separated continent:
	euler_pole_latitudes <- runif(length(separate_continents), -90, 90)
	
	# Ensure no latitude is directly at a pole (North or South) by redrawing if any are found:
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
		continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
		# If a negative or zero speed is picked then redraw:
		while(continent_speed <= 0) continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
		# Set degree change per step (effectively the speed):
		degrees_per_step[i] <- continent_speed / (2 * pi * furthest_continent_GC_distance) * 360
		
	}

	#starting information for the tree
	begin_cont <- ceiling(runif(1, 0, N_continents))
	life_begins <- EndPoint(continent_starting_points[begin_cont, "Longitude"], continent_starting_points[begin_cont, "Latitude"], runif(1,0,360), runif(1,0,radius))
	extra.rows <- matrix(NA,nrow=2,ncol= (N_steps * organism_multiplier) + 1)
	organism_lat_matrix <- matrix(nrow=2,ncol=(N_steps * organism_multiplier) + 1)
    organism_long_matrix <- matrix(nrow=2,ncol=(N_steps * organism_multiplier) + 1)
    organism_lat_matrix[,1] <- life_begins$lat
    organism_long_matrix[,1] <- life_begins$lon
    rownames(organism_lat_matrix) <- rep(begin_cont, 2)
    rownames(organism_long_matrix) <- rep(begin_cont, 2)
    edge <- rbind(c(1, 2), c(1, 3)) # this is a starting edge matrix
    edge.length <- rep(NA, 2)
    stem.depth <- numeric(2)
    alive <- rep(TRUE, 2) # marker for live lineages
    t <- 0; # time(step) at any point in the tree
    next.node <- 4

	# Now need to move them!
	#List to store which circles are in each supercontinent (new element added only when it changes)
	linked <- list()
	linked[[1]] <- separate_continents
	names(linked)[[1]] <- "1:1"
	
	#List to store which circles are in each supercontinent (new element added only when it changes)
	touching <- list()
	touching[[1]] <- touching_continents

	#Array to store the positions of every continent at each time step
	position <- array(NA, c(N_continents, N_steps + 1, 2), c("continent", "timestep", "coordinate"))
	position[,1,1] <- continent_starting_points[,"Longitude"]
	position[,1,2] <- continent_starting_points[,"Latitude"]

	temp_position<-matrix(nrow = N_continents, ncol=2)

	#for loops to move everything
	for (t in 2:(N_steps + 1)) {
		
		cat(t - 1, " ")
		
		#distances apart before they move
		starting_distances <- GreatCircleDistanceMatrix(position[,t-1,1], position[,t-1,2])

		for (k in 1:N_continents) {
			#find current longlat of the circle k
			start_long <- position[k, t-1, 1]
			start_lat <- position[k, t-1, 2]

			#Identify which supercontinent, and therefore which element of the euler pole and speed vectors, the circle k belongs to
			where <- which_supercontinent(k, tail(linked, n=1)[[1]])

			#Find distance of circle from pole
			distance <- GreatCircleDistanceFromLongLat(long1=start_long,lat1=start_lat, long2=euler_pole_longitudes[where], lat2=euler_pole_latitudes[where], Warn = FALSE)
			
			#Find bearing of circle from pole
			init_bearing <- BearingBetweenTwoLongLatPoints(euler_pole_longitudes[where], euler_pole_latitudes[where], start_long, start_lat)

			#Find the new bearing of circle from pole according to the speed specified
			new_bearing <- (init_bearing + degrees_per_step[where]) %% 360

			#Find the new location of the circle
			new_loc <- EndPoint(euler_pole_longitudes[where], euler_pole_latitudes[where], new_bearing, distance)

			#Add the new loction to the position matrix
			temp_position[k,1] <- new_loc$long
			temp_position[k,2] <- new_loc$lat

		}

		#Matrix of the distances between each continent in their new positions, before these are set
		new_distances <- GreatCircleDistanceMatrix(temp_position[,1], temp_position[,2])
		
		comp1 <- vector()
		comp2 <- vector()

		#Making a vector of whether any continents have collided and should be linked
		collisions<-matrix(nrow=0, ncol=2)
		for (b in 1:(N_continents-1)) {
			for (p in (b+1):N_continents) {
				comp1 <- all.equal(new_distances[p,b],starting_distances[p,b])
				comp2 <- new_distances[p,b] <= min_separation
				if (comp1 != TRUE && comp2 == TRUE) {
					collisions<-rbind(collisions,c(p,b))
				}
			}
		}
		
		perm_collisions <- matrix(nrow=0, ncol=2)

		# Moving continents back if there has only been one collision
		while(nrow(collisions) > 0) {

			# Set up vector to store proportional changes after collisions
			proportion <- vector()

			# Find out proportions to reduce to for all potential collisions
			for (coll in 1:nrow(collisions)) {
				cont_1 <- collisions[coll,1]
				cont_2 <- collisions[coll,2]
				where_1 <- which_supercontinent(cont_1, tail(linked, n=1)[[1]])
				where_2 <- which_supercontinent(cont_2, tail(linked, n=1)[[1]])
				continent_1_euler_longitude = euler_pole_longitudes[where_1]
				continent_1_euler_latitude = euler_pole_latitudes[where_1]
				continent_2_euler_longitude = euler_pole_longitudes[where_2]
				continent_2_euler_latitude = euler_pole_latitudes[where_2]
				continent_1_degrees_per_step = degrees_per_step[where_1]
				continent_2_degrees_per_step = degrees_per_step[where_2]
				
				proportion[coll] <- ColliderReverser(min_separation, position[cont_1, t-1, 1], position[cont_1, t-1, 2], temp_position[cont_1,1], temp_position[cont_1,2], position[cont_2, t-1, 1], position[cont_2, t-1, 2], temp_position[cont_2,1], temp_position[cont_2,2], continent_1_euler_longitude, continent_1_euler_latitude, continent_2_euler_longitude, continent_2_euler_latitude, continent_1_degrees_per_step, continent_2_degrees_per_step, EarthRad = 6367.4447, Warn = FALSE)
			
			}

			# Select proportion to move for first collision
			first_collision <- match(min(proportion),proportion)

			# Find the two continents that collided first
			cont_involved <- sort(collisions[first_collision, ])

			# Add to matrix of definite collisions that cannot be separated in the next step
			perm_collisions <- rbind(perm_collisions, cont_involved)

			# Move first clump back
			head_of_collision_1 <- cont_involved[1]
			where_1 <- which_supercontinent(head_of_collision_1, tail(linked, n=1)[[1]])

			#Find the other continents that were attached to the collider
			all_involved <- as.numeric(strsplit(tail(linked, n=1)[[1]][where_1], "&")[[1]])

			#Update temporary positions based on new change in bearing for all the continents in one of the clumps
			for (rev in 1:length(all_involved)) {
				cont_to_rev <- all_involved[rev]
				start_long <- position[cont_to_rev, t-1, 1]
				start_lat <- position[cont_to_rev, t-1, 2]
				distance <- GreatCircleDistanceFromLongLat(long1=start_long,lat1=start_lat, long2=euler_pole_longitudes[where_1], lat2=euler_pole_latitudes[where_1], Warn = FALSE)
			
				#Find bearing of circle from pole
				init_bearing <- BearingBetweenTwoLongLatPoints(euler_pole_longitudes[where_1], euler_pole_latitudes[where_1], start_long, start_lat)

				#Find the new bearing of circle from pole according to the speed specified
				new_bearing <- (init_bearing + (proportion[first_collision] * degrees_per_step[where_1])) %% 360

				#Find the new location of the circle
				new_loc <- EndPoint(euler_pole_longitudes[where_1], euler_pole_latitudes[where_1], new_bearing, distance)

				#Add the new loction to the position matrix
				temp_position[cont_to_rev,1] <- new_loc$long
				temp_position[cont_to_rev,2] <- new_loc$lat
				
			}

			#Move the second clump back
			head_of_collision_2 <- cont_involved[2]
			where_2 <- which_supercontinent(head_of_collision_2, tail(linked, n=1)[[1]])
			all_involved_2 <- as.numeric(strsplit(tail(linked, n=1)[[1]][where_2], "&")[[1]])

			#Update temporary positions based on new change in bearing for all the continents the other clump
			for (rev in 1:length(all_involved_2)) {
				cont_to_rev <- all_involved_2[rev]
				start_long <- position[cont_to_rev, t-1, 1]
				start_lat <- position[cont_to_rev, t-1, 2]
				distance <- GreatCircleDistanceFromLongLat(long1=start_long,lat1=start_lat, long2=euler_pole_longitudes[where_2], lat2=euler_pole_latitudes[where_2], Warn = FALSE)
			
				#Find bearing of circle from pole
				init_bearing <- BearingBetweenTwoLongLatPoints(euler_pole_longitudes[where_2], euler_pole_latitudes[where_2], start_long, start_lat)

				#Find the new bearing of circle from pole according to the speed specified
				new_bearing <-  (init_bearing + (proportion[first_collision] * degrees_per_step[where_2])) %% 360

				#Find the new location of the circle
				new_loc <- EndPoint(euler_pole_longitudes[where_2], euler_pole_latitudes[where_2], new_bearing, distance)

				#Add the new loction to the position matrix
				temp_position[cont_to_rev,1] <- new_loc$long
				temp_position[cont_to_rev,2] <- new_loc$lat
				
			}

			#Recalculate the new distances
			new_distances <- GreatCircleDistanceMatrix(temp_position[,1], temp_position[,2])
		
			comp1 <- vector()
			comp2 <- vector()

			# Check whether any other collisions have still occurred
			collisions <- matrix(nrow=0, ncol=2)
			
			for (b in 1:(N_continents - 1)) {
				
				for (p in (b+1):N_continents) {
					
					comp1 <- all.equal(new_distances[p, b], starting_distances[p, b])
					
					comp2 <- new_distances[p, b] <= min_separation
					
					if (comp1 != TRUE && comp2 == TRUE) {
					
						# If collision has not already been recorded:
						if(length(sort(match(paste(sort(c(b, p)), collapse=""), apply(perm_collisions, 1, paste, collapse="")))) == 0) {
						
							# Record collision:
							collisions <- rbind(collisions, c(p, b))
							
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
		separate_continents <- HowManySeparateContinents(min_separation, position[,t,1], position[,t,2])
		
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
						if(all(intersect(as.character(perm_collisions[j, ]), cluster_continent_numbers) == as.character(perm_collisions[j, ]))) {
						
							# Add new collision to protected links list:
							cluster_protected_links <- rbind(cluster_protected_links, perm_collisions[j, ])
							
						}
						
						
					}

				}
				
				# Get splits (new clusters) of separated cluster:
				splits <- ContinentSplitter(min_separation, cluster_longitudes, cluster_latitudes, cluster_continent_numbers, cluster_protected_links, EarthRad)
				
				# Update separate continents vector:
				separate_continents <- sort(c(separate_continents[-match(clusters_to_split[i], separate_continents)], splits))
				
			}
			
		}

# Now change the Euler poles and speeds
		
		# Get list of touching continents (to be used later for whether dispersal is allowable or not):
		touching_continents <- HowManySeparateContinents((radius * 2), position[,t,1], position[,t,2])
		if (any(touching_continents != tail(touching,n=1)[[1]])) {
			touching <- c(touching, list(touching_continents))
		}

		# If the continental clustering has changed:
		if (paste(sort(separate_continents), collapse="") != paste(sort(tail(linked,n=1)[[1]]), collapse="")) {
			
			# Add new continental configuration to linked list
			linked <- c(linked, list(separate_continents))

			# Update time step for previous continental configuration:
			names(linked)[(length(linked) - 1)] <- paste(strsplit(names(linked[(length(linked) - 1)]), ":")[[1]][1], ":", t - 1, sep="")
			
			# Add time step to new continental configuration:
			names(linked)[length(linked)] <- paste(t, ":", t, sep="")
			
			# Select continents that are different to previous time step
			toKeep <- c(tail(linked,n=1)[[1]], tail(linked, n=2)[[1]])[duplicated(c(tail(linked,n=1)[[1]], tail(linked, n=2)[[1]]))]

			#Vectors for new poles and speeds
			new_euler_latitudes <- rep(NA, length(separate_continents))
			new_euler_longitudes <- rep(NA, length(separate_continents))
			new_degrees_per_step <- rep(NA, length(separate_continents))

			if (length(toKeep) != 0){
				#Inherit previous poles and speeds for those clumps that remain the same
				for (y in 1:length(toKeep)) {
					new_euler_longitudes[match(toKeep[y],separate_continents)]<- euler_pole_longitudes[match(toKeep[y],tail(linked, n=2)[[1]])]
					new_euler_latitudes[match(toKeep[y],separate_continents)]<- euler_pole_latitudes[match(toKeep[y],tail(linked, n=2)[[1]])]
					new_degrees_per_step[match(toKeep[y],separate_continents)]<- degrees_per_step[match(toKeep[y],tail(linked, n=2)[[1]])]
				}
			}

			# Make new Euler poles
			new_euler_longitudes[is.na(new_euler_longitudes)]<- runif((length(separate_continents)-length(toKeep)), -180, 180)

			changed_euler_latitudes<- runif((length(separate_continents)-length(toKeep)), -90, 90)
			while((sum(changed_euler_latitudes == 90) + sum(changed_euler_latitudes == -90)) > 0) changed_euler_latitudes <- runif((length(separate_continents)-length(tokeep)), -90, 90)
			
			new_euler_latitudes[is.na(new_euler_latitudes)] <- changed_euler_latitudes

			# Get Great Circle distances from Euler pole to each continent centre:
			for (l in 1:(length(separate_continents)-length(toKeep))) {
				
				# Find the first NA to fill in:
				changer <- match(NA, new_degrees_per_step)
				
				#Find the contients that are in the clump whose speed is going to change
				rows <- as.numeric(unlist(strsplit(separate_continents[changer], "&")))
				
				#Find the new euler pole for that clump
				euler_long <- new_euler_longitudes[changer]
				euler_lat <- new_euler_latitudes[changer]
				
				#Find the distances of each of the continents in the clump from their new euler pole
				euler_GC_distances <- One2ManyGreatCircleDistance(euler_long, euler_lat, position[rows, t, 1], position[rows, t, 2])
				
				#Find which of those is closest to the equator
				furthest_continent_GC_distance <- (0.5 * pi * EarthRad) - min(abs(euler_GC_distances - (0.5 * pi * EarthRad)))
	
				# Randomly draw a continent speed:
				continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
				# If a negative or zero speed is picked then redraw:
				while(continent_speed <= 0) continent_speed <- rnorm(1, mean = continent_speed_mean, sd = continent_speed_sd)
		
				# Set degree change per step (effectively the speed):
				new_degrees_per_step[changer] <- continent_speed / (2 * pi * furthest_continent_GC_distance) * 360
			}
			organism_mover <- cbind()
			euler_pole_longitudes <- new_euler_longitudes
			euler_pole_latitudes <- new_euler_latitudes
			degrees_per_step <- new_degrees_per_step
		}
<<<<<<< HEAD

		#Move the animals with the continents
		for (cont in 1:N_continents) {
			moving <- which(rownames(organism_long_matrix)==cont)
			#Distanace the continent they're on moved
			for 
			
			#Find bearing of circle from pole
			init_bearing <- BearingBetweenTwoLongLatPoints(euler_pole_longitudes[where_2], euler_pole_latitudes[where_2], start_long, start_lat)
		}
=======
		
>>>>>>> origin/master
	}


	
	# Add final time step to continental configurations:
	names(linked)[length(linked)] <- paste(strsplit(names(linked)[length(linked)], ":")[[1]][1], ":", t, sep="")
	
	output <- list(position, linked)
	
	names(output) <- c("continent_positions", "continent_clusters")
	
	return(output)
	
}
# When rotating around Euler pole could theoretically pick clockwise or anticlockwise, but as we are allowing poles to be on either side of planet this takes care of that for us!
# Number continents in plots
	
#plot(position[1, , 1], position[1, , 2], xlim=c(-180, 180), ylim=c(-90, 90), col=rainbow(N_steps))

#for(i in 2:N_continents) points(position[i,,1], position[i,,2], col=rainbow(N_steps))
