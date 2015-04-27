#' Point(s) at which two Great Circles intersect on a sphere
#'
#' Given two Great Circles returns point(s) at which they intersect on a sphere.
#'
#' @param longitude_1 Decimalised longitude of first point on first Great Circle.
#' @param latitude_1 Decimalised latitude of first point on first Great Circle.
#' @param longitude_2 Decimalised longitude of second point on first Great Circle.
#' @param latitude_2 Decimalised latitude of second point on first Great Circle.
#' @param longitude_3 Decimalised longitude of first point on second Great Circle.
#' @param latitude_3 Decimalised latitude of first point on second Great Circle.
#' @param longitude_4 Decimalised longitude of second point on second Great Circle.
#' @param latitude_4 Decimalised latitude of second point on second Great Circle.
#' @param EarthRad Earth radius in kilometres.
#' @return Matrix of longitude-latitude points at which intersection(s) occur.
#' @details Nothing yet.
#' @examples
#' longitude_1 <- runif(1, -180, 180)
#' longitude_2 <- runif(1, -180, 180)
#' longitude_3 <- runif(1, -180, 180)
#' longitude_4 <- runif(1, -180, 180)
#' latitude_1 <- runif(1, -90, 90)
#' latitude_2 <- runif(1, -90, 90)
#' latitude_3 <- runif(1, -90, 90)
#' latitude_4 <- runif(1, -90, 90)
#' GreatCircleIntersection(longitude_1, latitude_1, longitude_2, latitude_2,
#'   longitude_3, latitude_3, longitude_4, latitude_4)

GreatCircleIntersection <- function(longitude_1, latitude_1, longitude_2, latitude_2, longitude_3, latitude_3, longitude_4, latitude_4, EarthRad = 6367.4447) {

# Set up input for passing to geosphere functions:
	
	p1 <- c(longitude_1, latitude_1)
	p2 <- c(longitude_2, latitude_2)
	p3 <- c(longitude_3, latitude_3)
	p4 <- c(longitude_4, latitude_4)
	
# All geosphere functions below this point:
	
# Author: Robert J. Hijmans
# April 2010
# version 1
# license GPL3
	
	.normalizeLonDeg <- function(x) {
		(x + 180) %% 360 - 180
	}
	
	.normalizeLonRad <- function(x) {
		(x + pi) %% (2*pi) - pi 
	}
	
	
	.isPolygon <- function(x, fix=FALSE) {
		x <- na.omit(x)
		if (nrow(x) < 4) {
			stop('this is not a polygon (insufficent number of vertices)')
		}
		if (length(unique(x[,1]))==1) {
			stop('All longitudes are the same (not a polygon)')
		}
		if (length(unique(x[,2]))==1) {
			stop('All latitudes are the same (not a polygon)')
		}
		if (! all(!(is.na(x))) ) {
			stop('polygon has NA values)')
		}
		if (! isTRUE(all.equal(x[1,], x[nrow(x),]))) {
			stop('this is not a valid (closed) polygon. The first vertex is not equal to the last vertex')	
		}
		return(x)
	}
	
# Author: Robert J. Hijmans
# October 2009
# version 1.0
# license GPL3
	
	antipodal <- function(p1, p2, tol=1e-9) {
		p1 <- .pointsToMatrix(p1) 
		p2 <- .pointsToMatrix(p2) 
		p <- cbind(p1[,1], p1[,2], p2[,1], p2[,2])	
		p[,c(1,3)] <- .normalizeLonDeg(p[,c(1,3)])
		diflon <- abs(p[,1] - p[,3]) 
		diflat <- abs(p[,2] + p[,4])
		(diflat < tol) & (diflon > (180 - tol))
	}
	
	
	antipode <- function(p) {
		p <- .pointsToMatrix(p)
		p[,1] <- .normalizeLonDeg(p[,1] + 180)
		p[,2] <- -p[,2]
		return( p )
	}
	
# Author: Robert J. Hijmans & Jacob van Etten
# October 2009
# version 1
# license GPL3
	
	.pointsToMatrix <- function(p, checkLonLat=TRUE, poly=FALSE) {
		if (inherits(p, 'SpatialPoints')) {
			test <- !is.projected(p)
			if (! isTRUE (test) ) {
				if (is.na(test)) {
					warning('Coordinate reference system of SpatialPoints object is not set. Assuming it is degrees (longitude/latitude)!')  			
				} else {
					stop('Points are projected. They should be in degrees (longitude/latitude)')  
				}
# or rather transform them ....?
			}
			p <- coordinates(p)
		} else if (is.data.frame(p)) {
			p <- as.matrix(p)
		} else 
		
		if (is.vector(p)){
			if (length(p) != 2) {
				stop('Wrong length for a vector, should be 2')
			} else {
				p <- matrix(p, ncol=2) 
			}
		} else if (is.matrix(p)) {
			if (length(p[1,]) != 2) {
				stop( 'A points matrix should have 2 columns')
			}
			cn <- colnames(p)
			if (length(cn) == 2) {
				if (toupper(cn[1]) == 'Y' | toupper(cn[2]) == 'X')  {
					warning('Suspect column names (x and y reversed?)')
				}
				if (toupper(substr(cn[1],1,3) == 'LAT' | toupper(substr(cn[2],1,3)) == 'LON'))  {
					warning('Suspect column names (longitude and latitude reversed?)')
				}
			}		
		} else {
			stop('points should be vectors of length 2, matrices with 2 columns, or inheriting from a SpatialPoints* object')
		}
		
		if (! is.numeric(p) ) { p[] <- as.numeric(p) }
		
		if (checkLonLat) {
			if (length(na.omit(p[,1])) > 0) {
				if (min(p[,1], na.rm=TRUE) < -360) { stop('longitude < -360') }
				if (max(p[,1], na.rm=TRUE) > 360) {  stop('longitude > 360')  }
				if (min(p[,1], na.rm=TRUE) < -180) { warning('longitude < -180') }
				if (max(p[,1], na.rm=TRUE) > 180) {  warning('longitude > 180')  }
			}
			if (length(na.omit(p[,2])) > 0) {
				if (min(p[,2], na.rm=TRUE) < -90) {  stop('latitude < -90')  }
				if (max(p[,2], na.rm=TRUE) > 90) {  stop('latitude > 90')  }
			}
		}
		
		
		if (poly) {
			if (! isTRUE(all.equal(p[1,], p[nrow(p),]))) {
				p <- rbind(p, p[1,])
			} 
			
			i <- p[-nrow(p),1] == p[-1,1] &  p[-nrow(p),2] == p[-1,2]
			i <- which(isTRUE(i))
			if (length(i) > 0) {
				p <- p[-i, ,drop=FALSE]
			}
			
			.isPolygon(p)
		}
		
		return(p)
	}
	
# author Robert Hijmans
# October 2009
# version 0.1
# license GPL3
	
# based on an alogrithm described by Ed Williams
# http://williams.best.vwh.net/intersect.htm
	
#intersection of two great circles defined by pt1 to pt2 and pt3 to pt4.
	
	einv <- function(e) {
		lat <- atan2(e[,3], sqrt(e[,1]^2 + e[,2]^2))
		lon <- atan2(-e[,2], e[,1]) 
		return(cbind(lon, lat))
	}
	
	eXe5 <- function(lon1, lat1, lon2, lat2) {
	    ex <- sin(lat1-lat2) *sin((lon1+lon2)/2) *cos((lon1-lon2)/2) - sin(lat1+lat2) *cos((lon1+lon2)/2) *sin((lon1-lon2)/2) 
		ey <- sin(lat1-lat2) *cos((lon1+lon2)/2) *cos((lon1-lon2)/2) + sin(lat1+lat2) *sin((lon1+lon2)/2) *sin((lon1-lon2)/2) 
		ez <- cos(lat1)*cos(lat2)*sin(lon1-lon2) 
		return( cbind(ex, ey, ez) )
	}
	
	eXe3 <- function(e1, e2) {
		x <- e1[,2] * e2[,3] -e2[,2] *e1[,3]
		y <- e1[,3] *e2[,1] -e2[,3] *e1[,1]
		z <- e1[,1] *e2[,2] -e1[,2] *e2[,1]
		return(cbind(x,y,z))
	}
	
	eSQRT <- function(e) {
		return(sqrt(e[,1]^2 + e[,2]^2 + e[,3]^2))
	}	
	
	p1 <- .pointsToMatrix(p1)
	p2 <- .pointsToMatrix(p2)
	p3 <- .pointsToMatrix(p3)
	p4 <- .pointsToMatrix(p4)
	
	p1 <- cbind(p1[,1], p1[,2], p2[,1], p2[,2])
	p3 <- cbind(p3[,1], p3[,2], p4[,1], p4[,2])
	p  <- cbind(p1[,1], p1[,2], p1[,3], p1[,4], p3[,1], p3[,2], p3[,3], p3[,4])
	
	p1 <- p[,1:2,drop=FALSE]
	p2 <- p[,3:4,drop=FALSE]
	p3 <- p[,5:6,drop=FALSE]
	p4 <- p[,7:8,drop=FALSE]
	
	res <- matrix(NA, nrow=nrow(p1), ncol=4)
	colnames(res) <- c('lon1', 'lat1', 'lon2', 'lat2')
	
	keep <- ! antipodal(p1, p2) | antipodal(p3, p4)
	keep <- keep & ! apply(p1 == p2, 1, sum) == 2
	
	if (sum(keep) == 0) { return(res) }
	
	toRad <- pi / 180 
	p1 <- p1[keep, , drop=FALSE] * toRad
	p2 <- p2[keep, , drop=FALSE] * toRad
	p3 <- p3[keep, , drop=FALSE] * toRad
	p4 <- p4[keep, , drop=FALSE] * toRad
	
	e1Xe2 <- eXe5(p1[,1], p1[,2], p2[,1], p2[,2])
	e3Xe4 <- eXe5(p3[,1], p3[,2], p4[,1], p4[,2])
	
	ea <- e1Xe2  / eSQRT(e1Xe2)
	eb <- e3Xe4  / eSQRT(e3Xe4)
	
	eaXeb <- eXe3(ea, eb)
	
	ll <- einv(eaXeb)
	ll2 <- cbind(ll[,1] + pi, -ll[,2])
	pts <- cbind(ll, ll2)
	pts[,1] <- .normalizeLonRad(pts[,1])
	pts[,3] <- .normalizeLonRad(pts[,3])
	
	res[keep,] <- pts / toRad

	return(res)
	
}
