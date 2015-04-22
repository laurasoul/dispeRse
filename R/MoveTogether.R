#' Move multiple 'continents' about a Euler pole
#'
#' Moves multiple 'continents' (which are temporarily or permanently joined) by their midpoints keeping them exactly the same distance from one another so that they behave as a single unit
#'
#' @param centers A matix with columns corresponding to continent central positions in degrees. Row 1 for longitude row 2 for latitude. Or a long lat vector for a single continent.
#' @param pole A vector of the long lat of the pole of rotation
#' @param angle The angle of rotation about the pole in degrees
#' @return a matix of the new long lat in degrees of the centers of each continent after moving
#' @keywords rotate sphere
#' @export
#' @examples
#' centers<-rbind(c(40, -50, -170), c(34, -60, 79))
#' pole<-c(90, -20)
#' MoveTogether(centers, pole, angle = 5)

# WHat happens when we try and break github?

MoveTogether <- function(centers, pole, angle) {
    
    if(pole[2] == 90 || pole[2] == -90) stop ("rotation pole cannot be at the north or south pole")
    
    new_centers<-matrix(nrow=2, ncol=ncol(centers))
    centers <- as.matrix(centers)
    
    for (i in 1:ncol(centers)) {

    	start_long <- centers[1,i]
    	start_lat <- centers[2,i]
    	distance <- GreatCircleDistanceFromLongLat(long1=start_long,lat1=start_lat, long2=pole[1], lat2=pole[2])
    	init_bearing <- BearingBetweenTwoLongLatPoints(pole[1], pole[2], start_long, start_lat)
    	new_bearing <- ((init_bearing + angle) + 360) %% 360
    	new_loc <- EndPoint(pole[1], pole[2], new_bearing, distance)
    	new_centers[1,i] <- new_loc$long
    	new_centers[2,i] <- new_loc$lat
    }

    return(new_centers)
}
