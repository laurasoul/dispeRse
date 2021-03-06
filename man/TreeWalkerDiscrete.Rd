% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TreeWalkerDiscrete.R
\name{TreeWalkerDiscrete}
\alias{TreeWalkerDiscrete}
\title{Generate random birth-death tree with associated coordinates}
\usage{
TreeWalkerDiscrete(
  b = 0.1,
  d = 0.05,
  steps = 50,
  slon = 0,
  slat = 0,
  steplengthsd = 100,
  EarthRad = 6367.4447
)
}
\arguments{
\item{b}{per-lineage birth (speciation) rate}

\item{d}{per-lineage death (extinction) rate}

\item{steps}{number of time steps to use}

\item{slon}{starting longitude}

\item{slat}{starting latitude}

\item{steplengthsd}{standard deviation used for random walk draws}

\item{EarthRad}{Earth radius in kilometres.}
}
\value{
tree a phylogenetic tree

longitudes a matrix with rows corresponding to the tree edges and colunns to time step

latitudes a matrix with rows corresponding to the tree edges and colunns to time step
}
\description{
This function generates a birth-death tree in discrete time steps at the same time as recording the long lat of each brach at the end of each step
}
\details{
This function is based on the function sim.bdtree in geiger <http://cran.r-project.org/web/packages/geiger/geiger.pdf>.
}
\examples{
TreeWalkerDiscrete(b=0.1, d=0.05, steps=50, slon=0, slat=0, steplengthsd = 100)
}
\author{
Laura C. Soul \email{lauracsoul@gmail.com}
}
\keyword{discrete}
\keyword{random}
\keyword{walk}
