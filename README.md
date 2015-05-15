dispeRse
========

This repository is an in development R package intended to [generate new models for detecting phylogenetic signals of endemism and dispersal](http://www.nescent.org/science/awards_summary.php?id=481) and is supported by two [Short-Term Visiting Scholarships](http://www.nescent.org/dir/visiting_scholar.php?type=Short-term%20Visitor) from [NESCent](http://www.nescent.org/) (to [Laura Soul](https://github.com/laurasoul) and [Graeme Lloyd](https://github.com/graemetlloyd)).

You can install and load dispeRse into R using the following:

```r
# Install the devtools package from CRAN:
install.packages("devtools")

# Load the devtools package into R:
library(devtools)

# Install the dispeRse package from github:
install_github("laurasoul/dispeRse")

# Load the dispeRse package into R:
library(dispeRse)
```

Script to perform (and plot) a random walk on a sphere:

```r
# Load libraries:
library(dispeRse)
install.packages("sphereplot", dependencies = TRUE)
library(sphereplot)

# Plot spherical grid:
rgl.sphgrid()

# Start at North Pole:
lonlat <- c(0, 90)

# For 1000 random walk steps:
for(i in 1:1000) {

	# Take a random step and store new coordinates:
	lonlat2 <- unlist(EndPoint(lonlat[1], lonlat[2], runif(1, 0, 360), abs(rnorm(1, 0, 100)))[c(1, 2)])

	# Plot coordinates on sphere in rainbow colour order:
	rgl.sphpoints(lonlat2[1], lonlat2[2], 1, deg=TRUE, col=rainbow(1000)[i], cex=2)

	# Update lonlat ready for next step:
	lonlat <- lonlat2

}
```

Script to perform (and plot) a random walk on a tree on a sphere in continuous time:

```r
# Add and load the Claddis package into R:
library(devtools)
install_github("graemetlloyd/Claddis")
library(Claddis)

# Load libraries:
library(dispeRse)
install.packages(c("ape", "maps"), dependencies = TRUE)
library(ape)
library(maps)

# Generate random 50-taxon tree:
tree <- rtree(10)

# Run function setting start point at equator-Greenwich Meridian intersection:
out <- TreeWalkerContinuous(tree, slon = 0, slat = 0, niter = 100, steplengthmean = 0, steplengthsd = 1000)

# Plot map:
map()

# For each branch:
for(i in 1:length(out)) {

	# Plot walk along branch:
	lines(out[[i]][2, ], out[[i]][3, ], col=rainbow(length(out))[i])

}
```

Script to perform (and plot) a random walk on a tree on a sphere in discrete time:

```r
# For discrete function
run <- TreeWalkerDiscrete(b = 0.1, d = 0.05, steps = 50, slon = 0, slat = 0, steplengthsd = 1000)

# Plot map:
map()

# For discrete function
for (i in 1:nrow(run$latitudes)) {

    lines(run$longitudes[i, ], run$latitudes[i, ], col = "red")

}
```

Script to check how different responses to hitting the edge of a continent affect eveness of distribution.

```r
# Load maps library:
library(maps)

# Set Earth radius in kilometres:
EarthRad <- 6367.4447

# Set number of time steps to run simulation for:
niter <- 100

# Set number of animals to run simulation for:
n_aminals <- 10000

# Set standard deviation of step size for dispersal in kilometres:
stepsize_sd <- 100

# Set centre of continent (long, lat):
continent_centre <- c(0, 0)

# Set radius of continent in kilometres:
continent_radius <- 2000

# Create empty matrix to store start positions (long, lat) of each animal:
start_positions <- matrix(nrow = 0, ncol = 2)

# Create vectors to store distances from continent centre at beginning and end of simlation:
end_distances <- start_distances <- vector(mode="numeric")

# Pick a response type (one of "bounce", "redraw", or "stop")
response <- "bounce"

# For each animal:
for(i in 1:n_aminals) {

	# Pick a random start position from the square in which the circular continent sits:
	start_position <- c(runif(1, -continent_radius / ((2 * pi * EarthRad) / 360), continent_radius / ((2 * pi * EarthRad) / 360)), runif(1, -continent_radius / ((2 * pi * EarthRad) / 360), continent_radius / ((2 * pi * EarthRad) / 360)))
	
	# Record distance from centre of continent:
	start_distances[i] <- GreatCircleDistanceFromLongLat(continent_centre[1], continent_centre[2], start_position[1], start_position[2])

	# If the starting point is not on the circular continent:
	if(start_distances[i] > continent_radius) {
	
		# Whilst the starting point is not on the circular continent:
		while(start_distances[i] > continent_radius) {
			
			# Pick a new starting point:
			start_position <- c(runif(1, -continent_radius / ((2 * pi * EarthRad) / 360), continent_radius / ((2 * pi * EarthRad) / 360)), runif(1, -continent_radius / ((2 * pi * EarthRad) / 360), continent_radius / ((2 * pi * EarthRad) / 360)))
			
			# Record distance from centre of continent:
			start_distances[i] <- GreatCircleDistanceFromLongLat(continent_centre[1], continent_centre[2], start_position[1], start_position[2])
			
		}
		
	}
	
	# Store point as start position (now we know it is safely on the continent:
	start_positions <- rbind(start_positions, start_position)
	
}

# Ossify start positions now as they will get overwritten later:
ossified_start_positions <- start_positions

# For each time step:
for(i in 1:niter) {
	
	# For each organism:
	for(j in 1:n_aminals) {
	
		# Draw a random bearing for the next move:
		bearing <- runif(1, 0, 360)
	
		# Draw a random distance for the next move:
		distance <- abs(rnorm(1, 0, stepsize_sd))
		
		# Calculate the new position of the organism following the move:
		new_position <- EndPoint(start_positions[j, 1], start_positions[j, 2], bearing, distance)[c("long", "lat")]

		# If draw means leaving continent:
		if(GreatCircleDistanceFromLongLat(new_position$long, new_position$lat, continent_centre[1], continent_centre[2]) > continent_radius) {
			
			# If response to leaving continent is to redraw:
			if(response == "redraw") {
				
				# As long as the draw removes the taxon from the continent:
				while(GreatCircleDistanceFromLongLat(new_position$long, new_position$lat, continent_centre[1], continent_centre[2]) > continent_radius) {
					
					# Redraw the bearing:
					bearing <- runif(1, 0, 360)
					
					# Redraw the distance:
					distance <- abs(rnorm(1, 0, stepsize_sd))
					
					# Calculate the new position:
					new_position <- EndPoint(start_positions[j, 1], start_positions[j, 2], bearing, distance)[c("long", "lat")]
					
				}
				
				# Overwrite start positions:
				start_positions[j, ] <- unlist(new_position)
				
				# Record distance from centre of continent:
				end_distances[j] <- GreatCircleDistanceFromLongLat(start_positions[j, 1], start_positions[j, 2], continent_centre[1], continent_centre[2])
				
			# If response to leaving continent is not to redraw (i.e., is to bounce off of, or stop at, continent edge):
			} else {
				
# FIND COLLISION POINT STARTS HERE:
				
				# Scalar which describes the proportion of the disatnce the organism has to travel before colliding with the continent edge:
				distance_modifier <- 0.5
				
				# Starting stepsize (will shrink as answer is honed in on more precisely):
				stepsize <- 0.1
				
				# Get starting new position (that will eb overwritten until the collision point is found):
				new_position <- EndPoint(start_positions[j, 1], start_positions[j, 2], bearing, distance_modifier * distance)[c("long", "lat")]
				
				# Whilst the collision point has not been found:
				while(all.equal(as.vector(GreatCircleDistanceFromLongLat(continent_centre[1], continent_centre[2], new_position$long, new_position$lat, EarthRad = EarthRad, Warn = FALSE)), continent_radius) != TRUE) {
					
					# Record current distance of new position from centre of continent:
					current_distance_from_centre <- as.vector(GreatCircleDistanceFromLongLat(continent_centre[1], continent_centre[2], new_position$long, new_position$lat, EarthRad = EarthRad, Warn = FALSE))
					
					# Do we need to increase the distance modifier value?:
					if(current_distance_from_centre < continent_radius) {
						
						# Create potential better value by adding step size to modifier:
						limit <- distance_modifier + stepsize
						
					# Or do we need to decrease the degree modifier value?:
					} else {
						
						# Create potential better value by subtracting step size from modifier:
						limit <- distance_modifier - stepsize
						
					}
					
					# Store the new new position (so we can ask if this is better):
					new_new_position <- EndPoint(start_positions[j, 1], start_positions[j, 2], bearing, limit * distance)[c("long", "lat")]
					
					# Get new distance from continent centre based on potential better modifer:
					new_distance_from_centre <- as.vector(GreatCircleDistanceFromLongLat(continent_centre[1], continent_centre[2], new_new_position$long, new_new_position$lat, EarthRad = EarthRad, Warn = FALSE))

					# If new distance is closer to continent radius:
					if(abs(current_distance_from_centre - continent_radius) > abs(new_distance_from_centre - continent_radius)) {
						
						# Update new position:
						new_position <- new_new_position
						
						# Update distance modifier itself:
						distance_modifier <- limit
						
					# If current distance is stil our best estimate:
					} else {
						
						# Shrink the step size so we can hone in closer:
						stepsize <- stepsize * 0.1
						
					}
					
				}
				
# FIND COLLISION POINT ENDS HERE:

				# If response to hitting a continent edge is to bounce off it:
				if(response == "bounce") {
					
					# Get the remaining distance of the step:
					remaining_distance <- as.vector(distance - GreatCircleDistanceFromLongLat(start_positions[j, 1], start_positions[j, 2], new_position$long, new_position$lat))
					
					# Get the bearing to the centre of the continent from the colision point:
					bearing_to_centre <- as.vector(BearingBetweenTwoLongLatPoints(new_position$long, new_position$lat, continent_centre[1], continent_centre[2]))
					
					# Get the bearing to the start of the dispersal step from the colision point:
					bearing_to_start <- as.vector(BearingBetweenTwoLongLatPoints(new_position$long, new_position$lat, start_positions[j, 1], start_positions[j, 2]))
					
					# Get the difference between these two bearings:
					bearing_difference <- min(c(abs(bearing_to_start - bearing_to_centre), 360 - abs(bearing_to_start - bearing_to_centre)))
					
					# If bearing difference should be subtracted to get new "bounce" bearing:
					if(all.equal((bearing_to_centre + bearing_difference) %% 360, bearing_to_start) == TRUE) {
						
						# Get new "bounce" bearing:
						bearing_to_end <- (bearing_to_centre - bearing_difference) %% 360
						
					# If bearing difference should be added to get new "bounce" bearing:
					} else {
						
						# Get new "bounce" bearing:
						bearing_to_end <- (bearing_to_centre + bearing_difference) %% 360
						
					}
					
					# Record actual new position (after bounce):
					new_position <- EndPoint(new_position$long, new_position$lat, bearing_to_end, remaining_distance)[c("long", "lat")]
					
					# Overwrite start positions:
					start_positions[j, ] <- unlist(new_position)
					
					# Record distance from centre of continent:
					end_distances[j] <- GreatCircleDistanceFromLongLat(start_positions[j, 1], start_positions[j, 2], continent_centre[1], continent_centre[2])
					
					# Little conditional to catch a second collision (does not do what it should, this is just a test to see if it happens):
					if(end_distances[j] > continent_radius) stop("Fallen off after bouncing!")
					
				}
				
				# If response to hitting a continent edge is to stick to it:
				if(response == "stop") {
					
					# Overwrite start positions:
					start_positions[j, ] <- unlist(new_position)
					
					# Record distance from centre of continent:
					end_distances[j] <- GreatCircleDistanceFromLongLat(start_positions[j, 1], start_positions[j, 2], continent_centre[1], continent_centre[2])
					
				}
				
			}
		
		# If draw means staying on continent:
		} else {
			
			# Overwrite start positions:
			start_positions[j, ] <- unlist(new_position)
			
			# Record distance from centre of continent:
			end_distances[j] <- GreatCircleDistanceFromLongLat(start_positions[j, 1], start_positions[j, 2], continent_centre[1], continent_centre[2])
			
		}
		
	}
	
}

# Now make plots to se if there is a bias in animal distribution based on response choice

# Plot map at a size just accommodating the circular continent:
map(xlim=c(-continent_radius / ((2 * pi * EarthRad) / 360), continent_radius / ((2 * pi * EarthRad) / 360)), ylim=c(-continent_radius / ((2 * pi * EarthRad) / 360), continent_radius / ((2 * pi * EarthRad) / 360)))

# Plot start positions in "red":
points(ossified_start_positions[, 1], ossified_start_positions[, 2], pch=19, col="red", cex=0.5)

# Plot end positions in blue:
points(start_positions[, 1], start_positions[, 2], pch=19, col="blue", cex=0.5)

# Create a set of 10 concentric rings to use in determining evenness of animal distribution:
concentric_circles <- seq(0, continent_radius, continent_radius / 10)

# Empty vectors to store N animals in each ring at beginning and end and the area of each ring:
start_N_animals_per_ring <- end_N_animals_per_ring <- sphere_areas <- vector(mode="numeric")

# For each ring:
for(i in 2:length(concentric_circles)) {
	
	# Get the area:
	sphere_areas[(i - 1)] <- SphericalCapArea(concentric_circles[i]) - SphericalCapArea(concentric_circles[(i - 1)])

	# Get the number of animals in the ring at the start of the simulation:
	start_N_animals_per_ring[(i - 1)] <- length(intersect(which(start_distances > concentric_circles[(i - 1)]), which(start_distances <= concentric_circles[i])))
	
	# Get the number of animals in the ring at the end of the simulation:
	end_N_animals_per_ring[(i - 1)] <- length(intersect(which(end_distances > concentric_circles[(i - 1)]), which(end_distances <= concentric_circles[i])))
	
}

# Make a barplot of N animals per area in each ring at the start of the simulation:
barplot(start_N_animals_per_ring / sphere_areas, main="Start")

# Make a barplot of N animals per area in each ring at the end of the simulation:
barplot(end_N_animals_per_ring / sphere_areas, main="End")
```

Animated results of a dispersal simulation run:

```r
# Load libraries (must have ImageMagick installed too):
library(dispeRse)
library(animation)

# Run a simulation (note, will crash if complete extinction occursso so re-run until it works):
out <- DispersalSimulator(N_steps = 1000, organism_multiplier = 1, radius = 2000, continent_speed_mean = 50, continent_speed_sd = 20, organism_step_sd = 1000, start_configuration = "supercontinent", stickiness = 0)

# Start outputting to GIF:
saveGIF({

	# For each time step:
	for(i in 1:1001) {
		
		# Plot continents:
		MapContinents(cbind(out$continent_positions[, i, 1], out$continent_positions[, i, 2]), radius = 2000)

		# Plot organisms:
		points(x = out$organism_longitudes[, i], y = out$organism_latitudes[, i], pch=19, col="red")

	}

# Additional output information:
}, movie.name = "test.gif", interval = 0.1, nmax = 1001, ani.width = 600, ani.height = 600)
```

More to follow.
