#' Which supercontinent?
#'
#' Which continental cluster does the current continent belong to.
#'
#' @param continent Continent number desired.
#' @param supercontinents Contienntla cluster vector.
#' @return Index of \code{continent} in \code{supercontinents}.
#' @details Nothing yet.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @export
#' @examples
#' WhichSupercontinent("3", c("1&2", "3&4"))

# Subfunction to find out which supercontinent a circle belongs to:
WhichSupercontinent <- function(continent, supercontinents) {
	
	# Make sure input data is in character format:
	continent <- as.character(continent)
	
	# Find out which continental cluster the continent belongs to:
	result <- which(unlist(lapply(lapply(lapply(strsplit(supercontinents, "&"), match, continent), sort), length)) == 1)
	
	# Output the result:
	return(result)
	
}
