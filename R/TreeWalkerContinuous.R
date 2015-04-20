#' Brownian motion (random walk) on a sphere
#'
#' Perform Brownian motion (a random walk) on a sphere.
#'
#' @param tree tree (as phylo object) to use
#' @param slon starting longitude
#' @param slat starting latitude
#' @param niter number of time steps to use
#' @param steplengthmean mean used for random walk draws
#' @param steplengthsd standard deviation used for random walk draws
#'
#' @return A list of edges with matrices showing times and coordinates
#' @keywords random walk
#' @export
#' @examples
#' tree <- rtree(10)
#' TreeWalkerContinuous(tree, slat = 0, slon = 0, niter = 1000, steplengthmean = 0, steplengthsd = 100)
#'

TreeWalkerContinuous <- function(tree, slon = 0, slat = 0, niter = 1000, steplengthmean = 0, steplengthsd = 100) {
	
# Add conditional to check tree is rooted!
# Add Claddis to dependencies?
# Hard to modify for changes in dispersal rate in future, not sure if can do anything about this though

	require(Claddis)
	
	# Establish root node number:
	root_node <- Ntip(tree) + 1
	
	# Get maximum path length and store as root age:
	tree$root.time <- max(diag(vcv(tree)))
	
	# Get step times:
	step_times <- seq(tree$root.time, 0, length.out=niter)

	# Get step size:
	step_size <- abs(diff(step_times)[1])
	
	# Get node ages:
	node_ages <- GetNodeAges(tree)
	
	# Get matrix of beginning and end ages for branches:
	edge_ages <- cbind(node_ages[tree$edge[, 1]], node_ages[tree$edge[, 2]], dimnames=c(c(), c()))
	
	# Matrix to store beginning and end points for each branch:
	end_points_matrix <- matrix(NA, ncol=4, nrow=nrow(tree$edge), dimnames=list(c(), c("Begin_lon", "Begin_lat", "End_lon", "End_lat")))
	
	# Get starting edges:	
	start_edges <- which(tree$edge[, 1] == root_node)
	
	# Create edges list:
	edges_list <- as.list(c(1:nrow(tree$edge)))
	
	# Sub-function to random walk along a single branch:
	BranchWalker <- function(slon = 0, slat = 0, step_times, step_size, edge_start_time, edge_end_time, steplengthmean = steplengthmean, steplengthsd = steplengthsd) {
	
		# Get first step present on branch:
		start_step <- min(which(edge_start_time >= step_times))
		
		# Get last step present on branch:
		end_step <- max(which(edge_end_time <= step_times))
		
		# Case if there is at least one step point along branch:
		if(start_step <= end_step) {
	
			# Case if branch starts at step value:
			if((edge_start_time - step_times[start_step]) == 0) {
				
				# Set starting longitude as input longitude:
				start_lon <- slon
				
				# Set starting longitude as input latitude:
				start_lat <- slat
				
			# Case if branch starts before a step value:
			} else {
				
				# Calculate fraction of a step being performed:
				step_fraction <- edge_start_time - step_times[start_step]
				
				# Make first step to join branch beginning to first step:
				first_step <- EndPoint(slong = slon, slat = slat, bearing = runif(1, 0, 360), distance = step_fraction * abs(rnorm(1, steplengthmean, steplengthsd)))
				
				# Record starting longitude:
				start_lon <- first_step$long
				
				# Record starting latitude:
				start_lat <- first_step$lat
				
			}
			
			# Create step matrix to store results:
			step_matrix <- matrix(nrow = 3, ncol = end_step - start_step + 1, dimnames=list(c("step_time", "step_long", "step_lat"), c(start_step:end_step)))
			
			# Fill in step times:
			step_matrix[1, ] <- step_times[start_step:end_step]
			
			# For each step:
			for(i in start_step:end_step) {
				
				# Record starting longitude:
				step_matrix[2, match(i, colnames(step_matrix))] <- start_lon
				
				# Record starting latitude:
				step_matrix[3, match(i, colnames(step_matrix))] <- start_lat
				
				# Make next step:
				next_step <- EndPoint(slong = start_lon, slat = start_lat, bearing = runif(1, 0, 360), distance = abs(rnorm(1, steplengthmean, steplengthsd)))
				
				# Update starting longitude:
				start_lon <- next_step$long
				
				# Update starting latitude:
				start_lat <- next_step$lat
				
			}
			
			# Case if branch ends at step value:
			if((step_times[end_step] - edge_end_time) == 0) {
				
				# Set ending longitude as last value in step matrix:
				elon <- step_matrix[2, as.character(end_step)]
				
				# Set ending latitude as last value in step matrix:
				elat <- step_matrix[3, as.character(end_step)]
				
			# Case if branch ends after a step value:
			} else {
				
				# Calculate fraction of a step from last step to end of branch:
				step_fraction <- (step_times[end_step] - edge_end_time) / step_size
				
				# Make actual last step to end of branch:
				last_step <- EndPoint(slong = step_matrix[2, as.character(end_step)], slat = step_matrix[3, as.character(end_step)], bearing = runif(1, 0, 360), distance = step_fraction * abs(rnorm(1, steplengthmean, steplengthsd)))
				
				# Record ending longitude:
				elon <- last_step$long
				
				# Record ending latitude:
				elat <- last_step$lat			
				
			}
			
			# Compile output:
			output <- cbind(c(edge_start_time, slon, slat), step_matrix, c(edge_end_time, elon, elat))
			
			# Update row names:
			rownames(output) <- c("step_time", "step_long", "step_lat")
			
			# Update column names:
			colnames(output)[1] <- "begin"
			
			# Update column names:
			colnames(output)[ncol(output)] <- "end"
			
		# Case if no steps on branch:
		} else {
			
			# Calculate fraction of a step represented by branch:
			step_fraction <- (edge_start_time - edge_end_time) / step_size
			
			# Make fractional step proportioned by branch length:
			branch_step <- EndPoint(slong = slon, slat = slat, bearing = runif(1, 0, 360), distance = step_fraction * abs(rnorm(1, steplengthmean, steplengthsd)))
			
			# Update starting longitude:
			elon <- branch_step$long
			
			# Update starting latitude:
			elat <- branch_step$lat
			
			# Compile output:
			output <- matrix(c(edge_start_time, edge_end_time, slon, elon, slat, elat), nrow=3, byrow=TRUE, dimnames=list(c("step_time", "step_long", "step_lat"), c("begin", "end")))
			
		}
			
		# Return output:
		return(output)
		
	}
	
	# For each starting (root connected) edge:
	for(i in start_edges) {
	
		# Get step_matrix for branch and add to edges_list:
		edges_list[[i]] <- BranchWalker(slon, slat, step_times, step_size, edge_ages[i, 1], edge_ages[i, 2], steplengthmean, steplengthsd)
		
		# Add end points of branches to end points matrix:
		end_points_matrix[i, ] <- c(edges_list[[i]][2:3, "begin"], edges_list[[i]][2:3, "end"])
		
	}

	# As long as there are branches that have not been walked along:
	while(any(is.na(end_points_matrix))) {
		
		# Find an edge that has starting point available:
		next_available_edge <- which(is.na(end_points_matrix[, 1]))[match(intersect(tree$edge[which(is.na(end_points_matrix[, 1])), 1], tree$edge[which(!is.na(end_points_matrix[, 1])), 2]), tree$edge[which(is.na(end_points_matrix[, 1])), 1])][1]
		
		# Find preceding edge (which has starting values for latitude and longitude at end):
		preceding_edge <- match(tree$edge[next_available_edge, 1], tree$edge[, 2])
		
		# Randomly walk along branch and store results:
		edges_list[[next_available_edge]] <- BranchWalker(edges_list[[preceding_edge]]["step_long", "end"], edges_list[[preceding_edge]]["step_lat", "end"], step_times, step_size, edge_ages[next_available_edge, 1], edge_ages[next_available_edge, 2], steplengthmean, steplengthsd)

		# Add end points of branches to end points matrix:
		end_points_matrix[next_available_edge, ] <- c(edges_list[[next_available_edge]][2:3, "begin"], edges_list[[next_available_edge]][2:3, "end"])

	}
	
	# Return edges list:
	return(edges_list)
	
}
