% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ThetaFromGreatCircleDistance.R
\name{ThetaFromGreatCircleDistance}
\alias{ThetaFromGreatCircleDistance}
\title{Get theta from great circl distance}
\usage{
ThetaFromGreatCircleDistance(great_circle_distance, EarthRad = 6367.4447)
}
\arguments{
\item{great_circle_distance}{The great circle distance in kilometres.}

\item{EarthRad}{Radius of the Earth in kilometres.}
}
\value{
Theta in radians.
}
\description{
Returns theta in radians from a great circle distance in kilometres
}
\details{
Nothing yet.
}
\examples{
ThetaFromGreatCircleDistance(1000)
}
\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}
