% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndPoint.R
\name{EndPoint}
\alias{EndPoint}
\title{Random walk on a sphere}
\usage{
EndPoint(
  longitude = 0,
  latitude = 0,
  bearing = 0,
  distance = 1,
  EarthRad = 6367.4447
)
}
\arguments{
\item{longitude}{starting longitude}

\item{latitude}{starting latitude}

\item{bearing}{heading in degrees from north}

\item{distance}{length of step in km}

\item{EarthRad}{Earth radius in kilometres.}
}
\value{
longitude final longitude

latitude final latitude

bearing final bearing in degrees from north
}
\description{
This function allows you calculate the final coordinates and bearing after one step of a random walk on a sphere.
}
\examples{
EndPoint(longitude = 0, latitude = 0, bearing = 90, distance = 111, EarthRad = 6367.4447)
}
\author{
Laura C. Soul \email{lauracsoul@gmail.com}
}
\keyword{random}
\keyword{walk}
