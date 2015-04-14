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
install.packages("sphereplot", dependencies=T)
library(sphereplot)

# Plot spherical grid:
rgl.sphgrid()

# Start at North Pole:
lonlat <- c(0, 90)

# For 1000 random walk steps:
for(i in 1:1000) {

	# Take a random step and store new coordinates:
	lonlat2 <- unlist(EndPoint(lonlat[2], lonlat[1], runif(1, 0, 360), abs(rnorm(1, 0, 100)))[c(2, 1)])

	# Plot coordinates on sphere in rainbow colour order:
	rgl.sphpoints(lonlat2[1], lonlat2[2], 1, deg=TRUE, col=rainbow(1000)[i], cex=2)

	# Update lonlat ready for next step:
	lonlat <- lonlat2

}
```

More to follow.
