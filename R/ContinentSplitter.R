#' Split apart joined continents
#'
#' 
#'
#' @param longitudes .
#' @param latitudes .
#' @return Continents.
#' @details Nothing yet.
#' @examples
#' 

ContinentSplitter <- function(min_separation, longitudes, latitudes, continent_numbers, protected_links, EarthRad = 6367.4447) {
	
	intercontinent_links <- IntercontinentalLinks(min_separation, longitudes, latitudes, EarthRad = EarthRad)
	
	
	
# Need to solve intersect problem first!
	
}
