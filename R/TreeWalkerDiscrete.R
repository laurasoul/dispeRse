#' Generate random birth-death tree with associated coordinates
#'
#' This function generates a birth-death tree in discrete time steps at the same time as recording the long lat of each brach at the end of each step
#'
#' @param b per-lineage birth (speciation) rate
#' @param d per-lineage death (extinction) rate
#' @param slon starting longitude
#' @param slat starting latitude
#' @param steps number of time steps to use
#' @param steplengthsd standard deviation used for random walk draws
#'
#' @return tree a phylogenetic tree
#' @return longitudes a matrix with rows corresponding to the tree edges and colunns to time step
#' @return latitudes a matrix with rows corresponding to the tree edges and colunns to time step
#' @details This function is based on the function sim.bdtree in geiger <http://cran.r-project.org/web/packages/geiger/geiger.pdf>. 
#' @keywords random walk discrete
#' @export
#' @examples
#' TreeWalkerDiscrete(b=0.1, d=0.05, steps=50, slon=0, slat=0, steplengthsd = 100)
#'

TreeWalkerDiscrete <- function (b=0.1, d=0.05, steps=50, slon=0, slat=0, steplengthsd=100) {
# Modified from sim.bdtree in geiger
# The following simulates birthdeath trees to a given number of time steps t
    extra.rows <- matrix(NA,nrow=2,ncol=steps+1)
    t <- steps
    time.stop = t
    if (time.stop == 0) stop("Stopping criterion ('n' or 't') must be provided")
    return.all.extinct = FALSE
    
    while (1) {
        lat.matrix <- matrix(nrow=2,ncol=steps+1)
        long.matrix <- matrix(nrow=2,ncol=steps+1)
        lat.matrix[,1] <- slat
        long.matrix[,1] <- slon
        edge <- rbind(c(1, 2), c(1, 3)) # this is a starting edge matrix
        edge.length <- rep(NA, 2)
        stem.depth <- numeric(2)
        alive <- rep(TRUE, 2) # marker for live lineages
        t <- 0; # time(step) at any point in the tree
        next.node <- 4
        
############
        repeat {
            if (sum(alive) == 0) break;
            dt<-1
            t <- t + dt
            for (i in 1:nrow(lat.matrix)) {
                if (alive[i]) {
                    starting<-c(lat.matrix[i,t],long.matrix[i,t])
                    moveto<-EndPoint(starting[2],starting[1],runif(1,0,360),abs(rnorm(1,0,steplengthsd))) #generates a random walk step and calculates new position
                    lat.matrix[i,t+1]<-moveto$lat
                    long.matrix[i,t+1]<-moveto$long
                }
            }

            if (time.stop) {
                if (t >= time.stop) {
                    t <- time.stop
                    break
                }
            }
            r <- runif(1)
            if (r <= b/(b + d)) { ###4 #this creates a bifucation in the tree
                random_lineage <- round(runif(1, min = 1, max = sum(alive)))
                e <- matrix(edge[alive,], ncol = 2)
                parent <- e[random_lineage,2]
                x <- which(edge[,2] == parent)
                new.rows.lat <- new.rows.long <- extra.rows
                new.rows.lat[, t+1] <- lat.matrix[x, t+1] #updates the long and lat matricies with new coordiantes for each lineage
                new.rows.long[, t+1] <- long.matrix[x, t+1]
                lat.matrix <- rbind(lat.matrix,new.rows.lat)
                long.matrix <- rbind(long.matrix,new.rows.long)
                alive[alive][random_lineage] <- FALSE
                edge <- rbind(edge, c(parent, next.node), c(parent, next.node + 1))
                next.node <- next.node + 2
                alive <- c(alive, TRUE, TRUE)
                stem.depth <- c(stem.depth, t, t)
                edge.length[x] <- t - stem.depth[x]
                edge.length<-c(edge.length, NA, NA)
            } else {###4 This terminates one of the current lineages on the tree
                random_lineage <- round(runif(1, min = 1, max = sum(alive)))
                edge.length[alive][random_lineage] <- t - stem.depth[alive][random_lineage]
                alive[alive][random_lineage] <- FALSE
            }###4
        }#1A
        
        if (return.all.extinct == TRUE | sum(alive) > 1) 
            break
    }
    edge.length[alive] <- t - stem.depth[alive];
    n <- -1
    for (i in 1:max(edge)) {
        if (any(edge[,1] == i)) {
            edge[which(edge[,1] == i), 1] <- n
            edge[which(edge[,2] == i), 2] <- n
            n <- n - 1;
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
    result <- list(obj,lat.matrix,long.matrix)
    names(result) <- c("tree", "latitudes", "longitudes")
    return (result)
}
